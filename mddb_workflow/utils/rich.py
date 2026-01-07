import contextlib

from rich.errors import NotRenderableError


# From: https://gist.github.com/izikeros/b0d32072f234fba73650eb4b1e9c0017
def rich_display_dataframe(df, title="Dataframe") -> None:
    """Display dataframe as table using rich library.

    Args:
        df (pd.DataFrame): dataframe to display
        title (str, optional): title of the table. Defaults to "Dataframe".

    Raises:
        NotRenderableError: if dataframe cannot be rendered

    Returns:
        rich.table.Table: rich table

    """
    from rich import print
    from rich.table import Table

    # ensure dataframe contains only string values
    df = df.astype(str)

    table = Table(title=title)
    for col in df.columns:
        table.add_column(col)
    for row in df.values:
        with contextlib.suppress(NotRenderableError):
            table.add_row(*row)
    print(table)
