from app import *
import logging

logger = logging.getLogger("panabyss_logger")

from dash import Input, Output, State, ctx, no_update, html, ALL
from sqlite_phylo_requests import *
import json


#Function used to read json params
def expand_jobs(df):

    # ---------
    # 1. parse params JSON
    # ---------
    def parse_params(x):
        try:
            return json.loads(x) if isinstance(x, str) else {}
        except:
            return {}

    params = df["params"].apply(parse_params).apply(pd.Series)

    # ---------
    # 2. flatten / sanitize list fields
    # ---------
    def safe_scalar(x):
        if isinstance(x, list):
            return ", ".join(map(str, x))
        if isinstance(x, dict):
            return json.dumps(x)
        return x

    df = df.drop(columns=["params"])
    df = pd.concat([df, params], axis=1)

    for col in df.columns:
        df[col] = df[col].apply(safe_scalar)

    return df

@app.callback(
    Output("phylo-jobs-table", "data"),
    Output("phylo-jobs-table", "columns"),
    Input("url", "pathname"),
)
def load_table(pathname):
    if pathname != "/phylogenetic_management":
        return no_update, no_update

    df = load_phylo_jobs()
    # 👉 expansion JSON + flatten lists
    df = expand_jobs(df)

    # Action column
    df["action"] = df["status"].apply(
        lambda s: "Cancel" if s == "RUNNING" else "Delete"
    )

    columns = [{"name": c, "id": c} for c in df.columns]

    return df.to_dict("records"), columns

@app.callback(
    Output("phylo-jobs-table", "data", allow_duplicate=True),
    Output("phylo-jobs-table", "columns", allow_duplicate=True),
    Input("phylo-jobs-table", "active_cell"),
    State("phylo-jobs-table", "data"),
    prevent_initial_call=True
)
def handle_action(active_cell, rows):

    if not active_cell or not rows:
        return no_update, no_update

    if active_cell["column_id"] != "action":
        return no_update, no_update

    job = rows[active_cell["row"]]
    job_id = job["job_id"]
    status = job["status"]

    if status == "RUNNING":
        set_phylo_job_cancel(job_id)
    else:
        delete_phylo_job(job_id)

    # reload
    df = load_phylo_jobs()
    df = expand_jobs(df)

    df["action"] = df["status"].apply(
        lambda s: "Cancel" if s == "RUNNING" else "Delete"
    )

    columns = [{"name": c, "id": c} for c in df.columns]

    return df.to_dict("records"), columns