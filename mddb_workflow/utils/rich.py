import contextlib
from rich.errors import NotRenderableError
from rich.live import Live
import time


# From: https://gist.github.com/izikeros/b0d32072f234fba73650eb4b1e9c0017
def rich_display_dataframe(df, row_limit=None, title="Dataframe", only_return=False) -> None:
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
        if row_limit and table.row_count >= row_limit:
            break
    if only_return:
        return table
    print(table)


def watch_dataframe(get_df_func, interval=2, title="Live Dataframe"):
    """Watch a dataframe with periodic updates."""
    with Live(refresh_per_second=0.5) as live:
        while True:
            df = get_df_func()
            table = rich_display_dataframe(
                df,
                row_limit=20,
                title=title,
                only_return=True
            )
            live.update(table)
            time.sleep(interval)
