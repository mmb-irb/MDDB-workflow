import contextlib
from rich.errors import NotRenderableError
from rich.live import Live
from rich.spinner import Spinner
from rich.console import Group
import time
import re


def _convert_html_links_to_rich(text: str) -> str:
    """Convert HTML links to Rich markup format.

    Args:
        text: String that may contain HTML anchor tags

    Returns:
        String with Rich markup links or original text

    """
    # Pattern to match <a href="URL" ...>text</a>
    pattern = r'<a href="([^"]+)"[^>]*>([^<]+)</a>'

    def replace_link(match):
        url = match.group(1)
        # link_text = match.group(2)
        return f'[link={url}]📝[/link]'

    result = re.sub(pattern, replace_link, text)
    return result


# From: https://gist.github.com/izikeros/b0d32072f234fba73650eb4b1e9c0017
def rich_display_dataframe(df, row_limit=None, title="Dataframe", only_return=False) -> None:
    """Display dataframe as table using rich library.

    Args:
        df (pd.DataFrame): dataframe to display
        row_limit (int, optional): maximum number of rows to display
        title (str, optional): title of the table. Defaults to "Dataframe".
        only_return (bool, optional): if True, return table without printing

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
        # Convert HTML links to Rich markup for each cell
        rich_row = [_convert_html_links_to_rich(cell) for cell in row]
        with contextlib.suppress(NotRenderableError):
            table.add_row(*rich_row)
        if row_limit and table.row_count >= row_limit:
            break
    if only_return:
        return table
    print(table)


def watch_dataframe(get_df_func, interval=2, title="Live Dataframe"):
    """Watch a dataframe with periodic updates."""
    spinner = Spinner("dots", text=title)
    with Live(refresh_per_second=20) as live:
        try:
            while True:
                df = get_df_func()
                table = rich_display_dataframe(
                    df,
                    row_limit=20,
                    title=spinner,
                    only_return=True
                )
                live.update(table)
                time.sleep(interval)
        except KeyboardInterrupt:
            pass
