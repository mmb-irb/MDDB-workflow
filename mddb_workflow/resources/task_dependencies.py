"""Build interactive dependency graphs for the MDDB workflow tasks.

The graph is derived from the live ``Task`` registries, task function
signatures, and declared output files.  Only relationships that cannot be
inferred from those registries are declared explicitly below: derived workflow
properties and project properties which forward values from a reference MD.
"""

from __future__ import annotations

from inspect import getfullargspec
from typing import Any, Literal

import networkx as nx


GraphScope = Literal['all', 'project', 'md']

GROUP_STYLES = {
    'Input files': '#2F855A',
    'Metadata inputs': '#38A169',
    'Configuration & flags': '#D69E2E',
    'Runtime services': '#718096',
    'Derived workflow data': '#319795',
    'Output files': '#DD6B20',
    'Preparation': '#3182CE',
    'Mapping & references': '#805AD5',
    'Analysis': '#E53E3E',
    'Metadata & outputs': '#D53F8C',
}

GROUP_ORDER = tuple(GROUP_STYLES)

ANALYSIS_FLAGS = {
    'apl',
    'channels',
    'clusters',
    'density',
    'dihedrals',
    'dist',
    'energies',
    'hbonds',
    'helical',
    'linter',
    'lorder',
    'markov',
    'pairwise',
    'pca',
    'perres',
    'pockets',
    'rgyr',
    'rmsds',
    'rmsf',
    'sas',
    'thickness',
    'tmscore',
}

MAPPING_FLAGS = {
    'chains',
    'inchikeys',
    'inchimap',
    'ligmap',
    'lipmap',
    'memmap',
    'pdbs',
    'protmap',
    'resmap',
}

OUTPUT_FLAGS = {'aiidata', 'mdmeta', 'pmeta', 'screenshot', 'stopology'}

INPUT_FILE_ARGUMENTS = {
    'aiida_data_file',
    'input_structure_file',
    'input_topology_file',
    'input_trajectory_files',
    'populations',
    'resorted_bonds_file',
    'resorted_charges_file',
    'transitions',
}

RUNTIME_SERVICE_ARGUMENTS = {
    'cache',
    'database',
    'parallel',
    'register',
    'remote',
    'warnings',
}

CONFIGURATION_ARGUMENTS = {
    'debug',
    'faith',
    'fit',
    'frames_limit',
    'guess_bonds',
    'ignore_bonds',
    'image',
    'interaction_cutoff',
    'interactions_auto',
    'local_blast',
    'maximum_pockets_number',
    'mercy',
    'most_populated_frames_number',
    'must_check_stable_bonds',
    'n_steps',
    'nodes_number',
    'overall_selection',
    'parameters',
    'patience',
    'projection_frames',
    'screenshot_frame',
    'time_splits',
    'translation',
    'verbose',
}

# These Project properties do not own Task descriptors.  They deliberately
# forward results produced for the reference MD (or for every MD in the case of
# ``md_charges``), so signature inspection alone cannot recover the link.
PROJECT_FORWARDED_MD_TASKS = {
    'charges': 'charges',
    'interactions': 'inter',
    'md_charges': 'charges',
    'snapshots': 'frames',
    'structure_file': 'inpro',
    'topology_file': 'inpro',
    'trajectory_file': 'inpro',
    'universe': 'mda_univ',
}

# These properties are calculated lazily rather than through ``Task``
# descriptors, so their relationships are not represented by task function
# signatures.  Keep the declarations close to the graph builder so the manual
# part of the model is explicit and easy to audit.
DERIVED_DEPENDENCIES = {
    'md': {
        'structure': ('structure_file',),
        'pbc_selection': ('structure', 'input_pbc_selection'),
        'pbc_residues': ('structure', 'pbc_selection'),
        'cg_selection': ('structure', 'input_cg_selection'),
        'cg_residues': ('structure', 'cg_selection'),
        'dummy_selection': ('structure', 'input_dummy_selection'),
        'forced_class_selections': (
            'structure',
            'input_forced_class_selections',
        ),
        'topology_reader': ('topology_file',),
        'dihedrals': ('topology_reader',),
    },
    'project': {
        'structure': ('structure_file',),
        'pbc_selection': ('structure', 'input_pbc_selection'),
        'pbc_residues': ('structure', 'pbc_selection'),
        'cg_selection': ('structure', 'input_cg_selection'),
        'cg_residues': ('structure', 'cg_selection'),
        'dummy_selection': ('structure', 'input_dummy_selection'),
        'topology_reader': ('topology_file',),
        'dihedrals': ('topology_reader',),
        'is_time_dependent': ('input_type',),
    },
}

COLLAPSIBLE_SOURCE_GROUPS = {
    'Configuration & flags',
    'Input files',
    'Metadata inputs',
    'Runtime services',
}


def _slug(value: str) -> str:
    return value.lower().replace(' & ', '-').replace(' ', '-').replace('_', '-')


def _pretty(value: str) -> str:
    return value.replace('_', ' ').strip().capitalize()


def _scope_label(scope: str) -> str:
    return 'MD' if scope == 'md' else 'Project'


def _task_node_id(scope: str, flag: str) -> str:
    return f'task:{scope}:{flag}'


def _task_group(flag: str) -> str:
    if flag in ANALYSIS_FLAGS:
        return 'Analysis'
    if flag in MAPPING_FLAGS:
        return 'Mapping & references'
    if flag in OUTPUT_FLAGS:
        return 'Metadata & outputs'
    return 'Preparation'


def _task_output_from_property_getter(
    getter: Any, task_type: type
) -> tuple[Any | None, str | None]:
    if isinstance(getter, task_type):
        return getter, None
    closure_values = {}
    free_variables = getattr(getattr(getter, '__code__', None), 'co_freevars', ())
    for name, cell in zip(free_variables, getattr(getter, '__closure__', ()) or ()):
        try:
            closure_values[name] = cell.cell_contents
        except ValueError:
            continue
    for value in closure_values.values():
        if isinstance(value, task_type):
            return value, closure_values.get('argument')
    return None, None


def _property_task_aliases(
    owner: type, task_type: type
) -> dict[str, dict[str, str | None]]:
    aliases = {}
    for property_name, value in vars(owner).items():
        if not isinstance(value, property) or value.fget is None:
            continue
        task, output_argument = _task_output_from_property_getter(
            value.fget, task_type
        )
        if task is not None:
            aliases[property_name] = {
                'flag': task.flag,
                'output_argument': output_argument,
            }
    return aliases


def _output_filename_label(output_filename: Any) -> str:
    if isinstance(output_filename, str):
        return output_filename
    return f'dynamic path: {getattr(output_filename, "__name__", repr(output_filename))}'


def _load_task_model() -> tuple[
    dict[str, dict[str, Any]],
    dict[str, dict[str, dict[str, str | None]]],
]:
    from mddb_workflow.mwf import MD, Project, md_requestables, project_requestables
    from mddb_workflow.utils.tasks import Task

    task_records: dict[str, dict[str, Any]] = {}
    registries = {'project': project_requestables, 'md': md_requestables}
    for scope, registry in registries.items():
        for flag, task in registry.items():
            if not isinstance(task, Task):
                continue
            try:
                arguments = tuple(getfullargspec(task.func).args)
            except TypeError:
                arguments = ()
            node_id = _task_node_id(scope, flag)
            task_records[node_id] = {
                'arguments': arguments,
                'fixed_arguments': frozenset(task.args),
                'flag': flag,
                'function': getattr(task.func, '__name__', repr(task.func)),
                'name': task.name,
                'output_arguments': frozenset(task.output_filenames),
                'output_filenames': {
                    argument: _output_filename_label(filename)
                    for argument, filename in task.output_filenames.items()
                },
                'scope': scope,
                'writes_output_directory': 'output_directory' in arguments,
            }

    aliases = {
        'project': _property_task_aliases(Project, Task),
        'md': _property_task_aliases(MD, Task),
    }
    for scope, property_aliases in aliases.items():
        for property_name, descriptor in property_aliases.items():
            node_id = _task_node_id(scope, str(descriptor['flag']))
            task_records[node_id].setdefault('output_properties', set()).add(
                property_name
            )
    for task in task_records.values():
        task['output_properties'] = tuple(sorted(task.get('output_properties', ())))
    return task_records, aliases


def _resolve_task_dependency(
    target_scope: str,
    argument: str,
    aliases: dict[str, dict[str, dict[str, str | None]]],
    task_records: dict[str, dict[str, Any]],
    output_nodes: dict[tuple[str, str, str], str],
) -> str | None:
    candidate_scopes = ('md', 'project') if target_scope == 'md' else ('project',)
    for candidate_scope in candidate_scopes:
        descriptor = aliases[candidate_scope].get(argument)
        flag = str(descriptor['flag']) if descriptor else None
        node_id = _task_node_id(candidate_scope, flag) if flag else None
        if node_id in task_records:
            output_argument = descriptor['output_argument'] if descriptor else None
            if output_argument:
                return output_nodes.get(
                    (candidate_scope, flag, output_argument), node_id
                )
            return node_id

    if target_scope == 'project':
        flag = PROJECT_FORWARDED_MD_TASKS.get(argument)
        node_id = _task_node_id('md', flag) if flag else None
        if node_id in task_records:
            descriptor = aliases['md'].get(argument)
            output_argument = descriptor['output_argument'] if descriptor else None
            if output_argument:
                return output_nodes.get(('md', flag, output_argument), node_id)
            return node_id
    return None


def _dependency_group(argument: str, task: dict[str, Any]) -> str:
    if argument in INPUT_FILE_ARGUMENTS:
        return 'Input files'
    if argument in RUNTIME_SERVICE_ARGUMENTS:
        return 'Runtime services'
    if argument.startswith('input_'):
        if task['flag'] in {'mdmeta', 'pmeta'}:
            return 'Metadata inputs'
        if argument.endswith(('_file', '_files', '_filepath', '_filepaths')):
            return 'Input files'
        return 'Configuration & flags'
    if argument in task['fixed_arguments'] or argument in CONFIGURATION_ARGUMENTS:
        return 'Configuration & flags'
    return 'Derived workflow data'


def _grouped_source_label(scope: str, group: str) -> str:
    suffixes = {
        'Configuration & flags': 'configuration & flags',
        'Input files': 'input files',
        'Metadata inputs': 'metadata inputs',
        'Runtime services': 'runtime services',
    }
    return f'{_scope_label(scope)} {suffixes[group]}'


def _add_dependency_edge(
    graph: nx.DiGraph,
    source: str,
    target: str,
    argument: str,
    relation: str = 'requires',
) -> None:
    if source == target:
        return
    if graph.has_edge(source, target):
        graph.edges[source, target]['arguments'].add(argument)
        graph.edges[source, target]['relations'].add(relation)
    else:
        graph.add_edge(
            source,
            target,
            arguments={argument},
            relations={relation},
        )


def _add_source_dependency(
    graph: nx.DiGraph,
    target: str,
    scope: str,
    argument: str,
    group: str,
    group_sources: bool,
) -> None:
    collapse = group_sources and group in COLLAPSIBLE_SOURCE_GROUPS
    if collapse:
        source = f'source:{scope}:{_slug(group)}'
        label = _grouped_source_label(scope, group)
        short_label = label
        kind = 'grouped dependency'
    else:
        source = f'source:{scope}:{argument}'
        label = _pretty(argument)
        short_label = argument
        kind = 'dependency'

    if source not in graph:
        graph.add_node(
            source,
            flag='',
            function='',
            group=group,
            kind=kind,
            label=label,
            members=set(),
            scope=scope,
            short_label=short_label,
        )
    graph.nodes[source]['members'].add(argument)
    _add_dependency_edge(graph, source, target, argument)


def _output_node_id(scope: str, flag: str, output_argument: str) -> str:
    return f'output:{scope}:{flag}:{output_argument}'


def _add_task_outputs(
    graph: nx.DiGraph,
    task_records: dict[str, dict[str, Any]],
    aliases: dict[str, dict[str, dict[str, str | None]]],
) -> dict[tuple[str, str, str], str]:
    output_nodes = {}
    for task_node, task in task_records.items():
        scope = task['scope']
        flag = task['flag']
        task_aliases = {
            property_name: descriptor
            for property_name, descriptor in aliases[scope].items()
            if descriptor['flag'] == flag
        }
        data_outputs = {
            property_name
            for property_name, descriptor in task_aliases.items()
            if descriptor['output_argument'] is None
        }
        graph.nodes[task_node]['declared_outputs'].update(data_outputs)

        output_arguments = set(task['output_arguments'])
        output_arguments.update(
            str(descriptor['output_argument'])
            for descriptor in task_aliases.values()
            if descriptor['output_argument'] is not None
        )
        for output_argument in sorted(output_arguments):
            property_names = sorted(
                property_name
                for property_name, descriptor in task_aliases.items()
                if descriptor['output_argument'] == output_argument
            )
            filename = task['output_filenames'].get(output_argument, '')
            static_filename = (
                filename if filename and not filename.startswith('dynamic path:') else ''
            )
            label = static_filename or (
                _pretty(property_names[0])
                if property_names
                else filename or _pretty(output_argument)
            )
            short_label = (
                static_filename
                or (property_names[0] if property_names else output_argument)
            )
            output_node = _output_node_id(scope, flag, output_argument)
            graph.add_node(
                output_node,
                filename=filename,
                flag='',
                function='',
                group='Output files',
                kind='output file',
                label=label,
                members=set(property_names),
                output_argument=output_argument,
                producer=task_node,
                scope=scope,
                short_label=short_label,
            )
            output_nodes[(scope, flag, output_argument)] = output_node
            declared_file_outputs = tuple(
                f'{property_name} — {filename}' if filename else property_name
                for property_name in property_names
            )
            graph.nodes[task_node]['declared_outputs'].update(
                declared_file_outputs or (filename or output_argument,)
            )
            _add_dependency_edge(
                graph,
                task_node,
                output_node,
                output_argument,
                relation='produces',
            )

        if task['writes_output_directory'] and not output_arguments:
            output_argument = 'output_directory'
            output_node = _output_node_id(scope, flag, output_argument)
            directory_label = f'{flag}/ output directory'
            graph.add_node(
                output_node,
                filename=f'{flag}/',
                flag='',
                function='',
                group='Output files',
                kind='output directory',
                label=directory_label,
                members=set(),
                output_argument=output_argument,
                producer=task_node,
                scope=scope,
                short_label=f'{flag}/',
            )
            output_nodes[(scope, flag, output_argument)] = output_node
            graph.nodes[task_node]['declared_outputs'].add(directory_label)
            _add_dependency_edge(
                graph,
                task_node,
                output_node,
                output_argument,
                relation='produces',
            )
    return output_nodes


def _derived_node_id(scope: str, property_name: str) -> str:
    return f'derived:{scope}:{property_name}'


def _resolve_known_dependency(
    scope: str,
    argument: str,
    aliases: dict[str, dict[str, dict[str, str | None]]],
    task_records: dict[str, dict[str, Any]],
    output_nodes: dict[tuple[str, str, str], str],
    derived_nodes: dict[tuple[str, str], str],
) -> str | None:
    derived_node = derived_nodes.get((scope, argument))
    if derived_node is not None:
        return derived_node
    return _resolve_task_dependency(
        scope,
        argument,
        aliases,
        task_records,
        output_nodes,
    )


def _add_derived_dependencies(
    graph: nx.DiGraph,
    task_records: dict[str, dict[str, Any]],
    aliases: dict[str, dict[str, dict[str, str | None]]],
    output_nodes: dict[tuple[str, str, str], str],
    group_sources: bool,
) -> dict[tuple[str, str], str]:
    derived_nodes = {
        (scope, property_name): _derived_node_id(scope, property_name)
        for scope, properties in DERIVED_DEPENDENCIES.items()
        for property_name in properties
    }
    for (scope, property_name), node_id in derived_nodes.items():
        graph.add_node(
            node_id,
            flag='',
            function='',
            group='Derived workflow data',
            kind='derived data',
            label=_pretty(property_name),
            members={property_name},
            scope=scope,
            short_label=property_name,
        )

    for scope, properties in DERIVED_DEPENDENCIES.items():
        for property_name, dependencies in properties.items():
            target = derived_nodes[(scope, property_name)]
            for argument in dependencies:
                source = _resolve_known_dependency(
                    scope,
                    argument,
                    aliases,
                    task_records,
                    output_nodes,
                    derived_nodes,
                )
                if source is not None:
                    _add_dependency_edge(
                        graph,
                        source,
                        target,
                        argument,
                        relation='derives',
                    )
                    continue
                group = _dependency_group(
                    argument,
                    {'fixed_arguments': frozenset(), 'flag': 'derived'},
                )
                _add_source_dependency(
                    graph,
                    target,
                    scope,
                    argument,
                    group,
                    group_sources,
                )
    return derived_nodes


def _freeze_attributes(graph: nx.DiGraph) -> None:
    for _, attributes in graph.nodes(data=True):
        attributes['members'] = tuple(sorted(attributes.get('members', ())))
        if 'declared_outputs' in attributes:
            attributes['declared_outputs'] = tuple(
                sorted(attributes['declared_outputs'])
            )
    for _, _, attributes in graph.edges(data=True):
        attributes['arguments'] = tuple(sorted(attributes.get('arguments', ())))
        attributes['relations'] = tuple(sorted(attributes.get('relations', ())))


def _remove_unconsumed_outputs(graph: nx.DiGraph) -> None:
    graph.remove_nodes_from(
        tuple(
            node
            for node, attributes in graph.nodes(data=True)
            if attributes['kind'] in {'output file', 'output directory'}
            and graph.out_degree(node) == 0
        )
    )


def _select_scope(graph: nx.DiGraph, scope: GraphScope) -> nx.DiGraph:
    if scope == 'all':
        return graph
    selected_tasks = {
        node
        for node, attributes in graph.nodes(data=True)
        if attributes['kind'] == 'task' and attributes['scope'] == scope
    }
    selected_nodes = set(selected_tasks)
    for task in selected_tasks:
        selected_nodes.update(nx.ancestors(graph, task))
        selected_nodes.update(
            target
            for _, target, edge in graph.out_edges(task, data=True)
            if 'produces' in edge['relations']
        )
    return graph.subgraph(selected_nodes).copy()


def _select_focus(graph: nx.DiGraph, focus: str | None) -> nx.DiGraph:
    if not focus:
        return graph
    matches = [
        node
        for node, attributes in graph.nodes(data=True)
        if attributes['kind'] in {'task', 'output file'}
        and (
            node == focus
            or attributes['kind'] == 'task'
            and attributes['flag'] == focus
            or attributes['kind'] == 'output file'
            and attributes['short_label'] == focus
        )
    ]
    if not matches:
        raise ValueError(f'Unknown task or output focus: {focus}')
    if len(matches) > 1:
        raise ValueError(
            f'Ambiguous task or output focus: {focus}. Use a full node id.'
        )
    selected = matches[0]
    selected_nodes = {selected, *nx.ancestors(graph, selected)}
    if graph.nodes[selected]['kind'] == 'task':
        selected_nodes.update(
            target
            for _, target, edge in graph.out_edges(selected, data=True)
            if 'produces' in edge['relations']
        )
    else:
        selected_nodes.update(nx.descendants(graph, selected))
    focused_graph = graph.subgraph(selected_nodes).copy()
    nx.set_node_attributes(focused_graph, False, 'focused')
    nx.set_node_attributes(focused_graph, False, 'directly_linked')
    focused_graph.nodes[selected]['focused'] = True
    directly_linked = {
        selected,
        *focused_graph.predecessors(selected),
        *focused_graph.successors(selected),
    }
    for node in directly_linked:
        focused_graph.nodes[node]['directly_linked'] = True
    nx.set_edge_attributes(
        focused_graph,
        {
            edge: selected in edge
            for edge in focused_graph.edges()
        },
        'directly_linked',
    )
    return focused_graph


def _refresh_group_members(graph: nx.DiGraph) -> None:
    for node, attributes in graph.nodes(data=True):
        if attributes['kind'] != 'grouped dependency':
            continue
        attributes['members'] = tuple(
            sorted(
                {
                    argument
                    for _, _, edge in graph.out_edges(node, data=True)
                    for argument in edge['arguments']
                }
            )
        )


def build_dependency_graph(
    scope: GraphScope = 'all',
    *,
    group_sources: bool = True,
    focus: str | None = None,
) -> nx.DiGraph:
    """Build a directed graph from workflow Task signatures.

    Edges point from a required task or data source to the task that consumes
    it.  Project and MD views retain cross-scope upstream tasks when they are
    real dependencies.  ``focus`` restricts the graph to one task and all of
    its ancestors.
    """
    if scope not in {'all', 'project', 'md'}:
        raise ValueError(f'Unsupported graph scope: {scope}')

    task_records, aliases = _load_task_model()
    graph = nx.DiGraph(scope=scope, grouped_sources=group_sources)
    ignored_arguments = {'output_directory', 'self', 'task'}
    for node_id, task in task_records.items():
        declared_inputs = tuple(
            argument
            for argument in task['arguments']
            if argument not in ignored_arguments
            and argument not in task['output_arguments']
        )
        graph.add_node(
            node_id,
            declared_inputs=declared_inputs,
            declared_outputs=set(),
            flag=task['flag'],
            function=task['function'],
            group=_task_group(task['flag']),
            kind='task',
            label=task['name'],
            members=(),
            scope=task['scope'],
            short_label=task['flag'],
        )

    output_nodes = _add_task_outputs(graph, task_records, aliases)
    derived_nodes = _add_derived_dependencies(
        graph,
        task_records,
        aliases,
        output_nodes,
        group_sources,
    )
    for target, task in task_records.items():
        for argument in task['arguments']:
            if argument in ignored_arguments or argument in task['output_arguments']:
                continue
            source_task = _resolve_known_dependency(
                task['scope'],
                argument,
                aliases,
                task_records,
                output_nodes,
                derived_nodes,
            )
            if source_task is not None:
                _add_dependency_edge(graph, source_task, target, argument)
                continue
            group = _dependency_group(argument, task)
            _add_source_dependency(
                graph,
                target,
                task['scope'],
                argument,
                group,
                group_sources,
            )

    _remove_unconsumed_outputs(graph)
    _freeze_attributes(graph)
    scoped_graph = _select_scope(graph, scope)
    selected_graph = _select_focus(scoped_graph, focus)
    _refresh_group_members(selected_graph)
    return selected_graph


def task_options(scope: GraphScope = 'all') -> dict[str, str]:
    """Return task and connected-output choices for focus controls."""
    if scope not in {'all', 'project', 'md'}:
        raise ValueError(f'Unsupported graph scope: {scope}')
    task_records, _ = _load_task_model()
    options = {}
    for node_id, task in sorted(
        task_records.items(), key=lambda item: (item[1]['scope'], item[1]['flag'])
    ):
        if scope != 'all' and task['scope'] != scope:
            continue
        label = f'{_scope_label(task["scope"])} · {task["flag"]} — {task["name"]}'
        options[label] = node_id

    graph = build_dependency_graph(scope)
    output_nodes = sorted(
        (
            node,
            attributes,
        )
        for node, attributes in graph.nodes(data=True)
        if attributes['kind'] == 'output file'
    )
    for node_id, attributes in output_nodes:
        producer = graph.nodes[attributes['producer']]['flag']
        label = (
            f'{_scope_label(attributes["scope"])} · output — '
            f'{attributes["short_label"]} (from {producer})'
        )
        options[label] = node_id
    return options


def graph_statistics(graph: nx.DiGraph) -> dict[str, Any]:
    """Return compact counts for display alongside a dependency graph."""
    task_nodes = [
        node
        for node, attributes in graph.nodes(data=True)
        if attributes['kind'] == 'task'
    ]
    output_nodes = [
        node
        for node, attributes in graph.nodes(data=True)
        if attributes['kind'] in {'output file', 'output directory'}
    ]
    task_groups = {
        group: sum(graph.nodes[node]['group'] == group for node in task_nodes)
        for group in GROUP_ORDER
        if any(graph.nodes[node]['group'] == group for node in task_nodes)
    }
    return {
        'dependencies': graph.number_of_nodes() - len(task_nodes) - len(output_nodes),
        'edges': graph.number_of_edges(),
        'outputs': len(output_nodes),
        'task_groups': task_groups,
        'tasks': len(task_nodes),
    }


def _layered_positions(graph: nx.DiGraph) -> tuple[dict[str, tuple[float, float]], int]:
    if not graph:
        return {}, 0

    condensed = nx.condensation(graph)
    component_layers: dict[int, int] = {}
    for component in nx.topological_sort(condensed):
        predecessors = list(condensed.predecessors(component))
        component_layers[component] = (
            max(component_layers[parent] for parent in predecessors) + 1
            if predecessors
            else 0
        )

    node_layers = {
        node: component_layers[component]
        for node, component in condensed.graph['mapping'].items()
    }
    layers: dict[int, list[str]] = {}
    for node, layer in node_layers.items():
        layers.setdefault(layer, []).append(node)

    positions = {}
    max_layer_size = 0
    group_rank = {group: index for index, group in enumerate(GROUP_ORDER)}
    for layer, nodes in layers.items():
        nodes.sort(
            key=lambda node: (
                group_rank[graph.nodes[node]['group']],
                graph.nodes[node]['scope'],
                graph.nodes[node]['short_label'],
            )
        )
        max_layer_size = max(max_layer_size, len(nodes))
        midpoint = (len(nodes) - 1) / 2
        for index, node in enumerate(nodes):
            positions[node] = (layer * 3.2, (midpoint - index) * 1.35)
    return positions, max_layer_size


def _hover_items(title: str, items: tuple[str, ...] | list[str]) -> str:
    values = tuple(items)
    lines = '<br>'.join(f'&nbsp;&nbsp;• {item}' for item in values)
    return f'<br><b>{title}</b><br>{lines or "&nbsp;&nbsp;None declared"}'


def _hover_text(graph: nx.DiGraph, node: str) -> str:
    attributes = graph.nodes[node]
    scope = _scope_label(attributes['scope'])
    if attributes['kind'] == 'task':
        return (
            f'<b>{attributes["flag"]} — {attributes["label"]}</b><br>'
            f'Group: {attributes["group"]}<br>'
            f'Level: {scope}<br>'
            f'Function: {attributes["function"]}'
            f'{_hover_items("Inputs", attributes["declared_inputs"])}'
            f'{_hover_items("Outputs", attributes["declared_outputs"])}'
        )

    inputs = sorted(
        graph.nodes[source]['short_label'] for source in graph.predecessors(node)
    )
    outputs = sorted(
        graph.nodes[target]['short_label'] for target in graph.successors(node)
    )
    members = attributes.get('members', ())
    details = ''
    if attributes['kind'] in {'output file', 'output directory'}:
        details = (
            f'<br>Filename: {attributes["filename"]}'
            f'<br>Output argument: {attributes["output_argument"]}'
        )
    elif members:
        details = _hover_items('Contains', members)
    return (
        f'<b>{attributes["label"]}</b><br>'
        f'Group: {attributes["group"]}<br>'
        f'Level: {scope}'
        f'{details}'
        f'{_hover_items("Inputs", inputs)}'
        f'{_hover_items("Outputs / used by", outputs)}'
    )


def create_dependency_figure(graph: nx.DiGraph, *, title: str | None = None) -> Any:
    """Create an interactive Plotly figure for a dependency graph."""
    try:
        import plotly.graph_objects as go
    except ImportError as error:
        raise RuntimeError(
            'Plotly is required for this visualization. Install it with '
            '`conda install -c conda-forge plotly` or `pip install plotly`.'
        ) from error

    positions, max_layer_size = _layered_positions(graph)
    figure = go.Figure()
    focus_active = any(
        attributes.get('focused') for _, attributes in graph.nodes(data=True)
    )
    focused_node = next(
        (
            node
            for node, attributes in graph.nodes(data=True)
            if attributes.get('focused')
        ),
        'all',
    )

    for group in GROUP_ORDER:
        group_edges = [
            (source, target)
            for source, target in graph.edges()
            if graph.nodes[source]['group'] == group
        ]
        edge_sets = (
            (
                (
                    direct,
                    [
                        edge
                        for edge in group_edges
                        if graph.edges[edge].get('directly_linked', False) is direct
                    ],
                )
                for direct in (False, True)
            )
            if focus_active
            else ((None, group_edges),)
        )
        for direct, visible_edges in edge_sets:
            if not visible_edges:
                continue
            edge_x: list[float | None] = []
            edge_y: list[float | None] = []
            arrow_x = []
            arrow_y = []
            for source, target in visible_edges:
                source_x, source_y = positions[source]
                target_x, target_y = positions[target]
                edge_x.extend((source_x, target_x, None))
                edge_y.extend((source_y, target_y, None))
                arrow_x.append(source_x + (target_x - source_x) * 0.84)
                arrow_y.append(source_y + (target_y - source_y) * 0.84)
            line_opacity = 1.0 if direct else 0.08 if focus_active else 0.42
            arrow_opacity = 1.0 if direct else 0.1 if focus_active else 0.78
            figure.add_trace(
                go.Scatter(
                    x=edge_x,
                    y=edge_y,
                    mode='lines',
                    line={'color': GROUP_STYLES[group], 'width': 1.1},
                    opacity=line_opacity,
                    hoverinfo='skip',
                    legendgroup=group,
                    showlegend=False,
                )
            )
            figure.add_trace(
                go.Scatter(
                    x=arrow_x,
                    y=arrow_y,
                    mode='markers',
                    marker={
                        'color': GROUP_STYLES[group],
                        'size': 7,
                        'symbol': 'triangle-right',
                    },
                    opacity=arrow_opacity,
                    hoverinfo='skip',
                    legendgroup=group,
                    showlegend=False,
                )
            )

        nodes = [
            node for node, attributes in graph.nodes(data=True) if attributes['group'] == group
        ]
        if not nodes:
            continue
        figure.add_trace(
            go.Scatter(
                x=[positions[node][0] for node in nodes],
                y=[positions[node][1] for node in nodes],
                mode='markers+text',
                name=group,
                legendgroup=group,
                text=[graph.nodes[node]['short_label'] for node in nodes],
                textposition='bottom center',
                textfont={
                    'color': [
                        'rgba(26, 32, 44, 1)'
                        if graph.nodes[node].get('directly_linked')
                        else 'rgba(26, 32, 44, 0.12)'
                        if focus_active
                        else '#1A202C'
                        for node in nodes
                    ],
                    'size': 11,
                },
                customdata=[_hover_text(graph, node) for node in nodes],
                hovertemplate='%{customdata}<extra></extra>',
                cliponaxis=False,
                marker={
                    'color': GROUP_STYLES[group],
                    'line': {
                        'color': [
                            '#1A202C' if graph.nodes[node].get('focused') else '#FFFFFF'
                            for node in nodes
                        ],
                        'width': [
                            4 if graph.nodes[node].get('focused') else 1.5 for node in nodes
                        ],
                    },
                    'opacity': [
                        1.0
                        if graph.nodes[node].get('directly_linked')
                        else 0.12
                        if focus_active
                        else 0.94
                        for node in nodes
                    ],
                    'size': [
                        34
                        if graph.nodes[node].get('focused')
                        else 25
                        if graph.nodes[node]['kind'] == 'task'
                        else 20
                        for node in nodes
                    ],
                    'symbol': [
                        'circle'
                        if graph.nodes[node]['kind'] == 'task'
                        else 'square'
                        if graph.nodes[node]['kind'] == 'grouped dependency'
                        else 'hexagon'
                        if graph.nodes[node]['kind']
                        in {'output file', 'output directory'}
                        else 'diamond'
                        for node in nodes
                    ],
                },
            )
        )

    scope = graph.graph.get('scope', 'all')
    default_titles = {
        'all': 'MDDB workflow task dependencies',
        'project': 'Project task dependencies',
        'md': 'MD task dependencies',
    }
    figure.update_layout(
        title={'text': title or default_titles[scope], 'x': 0.02, 'xanchor': 'left'},
        template='plotly_white',
        height=max(680, min(1500, max_layer_size * 64 + 220)),
        margin={'b': 45, 'l': 35, 'r': 35, 't': 130},
        dragmode='pan',
        hovermode='closest',
        legend={
            'groupclick': 'togglegroup',
            'orientation': 'h',
            'title': {'text': 'Dependency groups'},
            'x': 0,
            'y': 1.12,
        },
        xaxis={
            'fixedrange': False,
            'showgrid': False,
            'showticklabels': False,
            'title': 'Dependencies  →  execution order',
            'zeroline': False,
        },
        yaxis={
            'fixedrange': False,
            'showgrid': False,
            'showticklabels': False,
            'zeroline': False,
        },
        # Reset Plotly's retained viewport when the picker changes the graph.
        uirevision=(
            f'{scope}-{graph.graph.get("grouped_sources", True)}-{focused_node}'
        ),
    )
    return figure


__all__ = [
    'GROUP_STYLES',
    'build_dependency_graph',
    'create_dependency_figure',
    'graph_statistics',
    'task_options',
]
