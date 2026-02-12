import os
import glob
import time
import traceback
import requests

# -----------------------------
# Config
# -----------------------------
CYREST = "http://127.0.0.1:1234/v1"

EDGE_DIR = os.path.join("07_network_construction", "edge_lists")
NODE_DIR = os.path.join("07_network_construction", "nodes_attributes")
PNG_DIR = os.path.join("07_network_construction", "cytoscape_pngs")
os.makedirs(PNG_DIR, exist_ok=True)

STYLE_NAME = "auto_group_style"
LAYOUT_NAME = "force-directed"  # fallback if not supported

FEATURE_COLOR = "#4C78A8"
CLUSTER_COLOR = "#F58518"
FEATURE_SHAPE = "ELLIPSE"
CLUSTER_SHAPE = "DIAMOND"

PNG_WIDTH = 2400
PNG_HEIGHT = 1800

# -----------------------------
# REST helpers
# -----------------------------
def _url(path: str) -> str:
    return CYREST + path

def cy_get(path: str):
    r = requests.get(_url(path))
    r.raise_for_status()
    return r.json() if r.text else None

def cy_post(path: str, **kwargs):
    r = requests.post(_url(path), **kwargs)
    r.raise_for_status()
    return r.json() if r.text else None

def cy_delete(path: str):
    r = requests.delete(_url(path))
    r.raise_for_status()
    return r.json() if r.text else None

def cytoscape_ping():
    try:
        cy_get("/version")
        return True
    except Exception:
        return False

def prefix_from_edge_path(edge_path: str) -> str:
    base = os.path.basename(edge_path)
    return base[:-len("_edge_list.csv")] if base.endswith("_edge_list.csv") else os.path.splitext(base)[0]

# -----------------------------
# Cytoscape actions
# -----------------------------
def import_network_from_csv(edge_csv: str, title: str) -> int:
    params = {
        "title": title,
        "sourceColumn": "node1",
        "targetColumn": "node2",
        "firstRowAsColumnNames": "true",
    }

    with open(edge_csv, "rb") as f:
        files = {"file": (os.path.basename(edge_csv), f, "text/csv")}
        res = cy_post("/networks/file", params=params, files=files)

    # Normalize response to an int network SUID
    if isinstance(res, int):
        net_suid = res
    elif isinstance(res, dict) and "networkSUID" in res:
        net_suid = int(res["networkSUID"])
    elif isinstance(res, list) and len(res) > 0:
        net_suid = int(res[0]) if not isinstance(res[0], dict) else int(res[0].get("networkSUID", res[0].get("suid")))
    else:
        raise RuntimeError(f"Unexpected network import response: {res}")

    # Ensure a view exists
    views = cy_get(f"/networks/{net_suid}/views")
    if not views:
        cy_post(f"/networks/{net_suid}/views")
        time.sleep(0.2)

    return net_suid

def import_node_table(node_csv: str, net_suid: int):
    # Map CSV col 0 (Node) to Cytoscape node key "name"
    params = {
        "keyColumnIndex": "0",
        "keyColumnForMapping": "name",
        "firstRowAsColumnNames": "true"
    }

    with open(node_csv, "rb") as f:
        files = {"file": (os.path.basename(node_csv), f, "text/csv")}
        cy_post(f"/networks/{net_suid}/tables/defaultnode", params=params, files=files)

def get_style_names() -> set:
    """
    Cytoscape can return:
      - ["default", "foo", ...] OR
      - [{"title":"default"}, {"title":"foo"}, ...]
    """
    styles = cy_get("/styles")
    if not styles:
        return set()

    if isinstance(styles, list):
        # list of strings
        if all(isinstance(x, str) for x in styles):
            return set(styles)
        # list of dicts
        if all(isinstance(x, dict) for x in styles):
            return {x.get("title") for x in styles if x.get("title")}
    # fallback: unknown format
    return set()

def ensure_style():
    if STYLE_NAME in get_style_names():
        return

    style_def = {
        "title": STYLE_NAME,
        "defaults": [
            {"visualProperty": "NODE_SIZE", "value": 30},
            {"visualProperty": "EDGE_WIDTH", "value": 2},
            {"visualProperty": "NODE_LABEL_FONT_SIZE", "value": 10},
        ],
        "mappings": [
            {
                "mappingType": "discrete",
                "mappingColumn": "Group",
                "mappingColumnType": "String",
                "visualProperty": "NODE_FILL_COLOR",
                "map": [
                    {"key": "feature", "value": FEATURE_COLOR},
                    {"key": "cluster", "value": CLUSTER_COLOR},
                ],
            },
            {
                "mappingType": "discrete",
                "mappingColumn": "Group",
                "mappingColumnType": "String",
                "visualProperty": "NODE_SHAPE",
                "map": [
                    {"key": "feature", "value": FEATURE_SHAPE},
                    {"key": "cluster", "value": CLUSTER_SHAPE},
                ],
            },
        ],
    }

    cy_post("/styles", json=style_def)

def apply_style_to_network(net_suid: int):
    # Some versions support this endpoint:
    cy_post(f"/apply/styles/{STYLE_NAME}/{net_suid}")

def layout_network(net_suid: int):
    # Try preferred layout then fallback
    try:
        cy_post(f"/apply/layouts/{LAYOUT_NAME}/{net_suid}")
    except Exception:
        cy_post(f"/apply/layouts/prefuse-force-directed/{net_suid}")

def export_png(net_suid: int, out_png: str):
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    r = requests.get(_url(f"/networks/{net_suid}/views/first.png"),
                     params={"h": PNG_HEIGHT, "w": PNG_WIDTH})
    r.raise_for_status()
    with open(out_png, "wb") as f:
        f.write(r.content)

def delete_network(net_suid: int):
    try:
        cy_delete(f"/networks/{net_suid}")
    except Exception:
        pass

# -----------------------------
# Main
# -----------------------------
def main():
    if not cytoscape_ping():
        raise RuntimeError(
            "Cannot reach Cytoscape at http://127.0.0.1:1234.\n"
            "1) Open Cytoscape\n"
            "2) Ensure CyREST is enabled (default)\n"
            "3) Check port in Cytoscape preferences if needed"
        )

    ensure_style()

    edge_files = sorted(glob.glob(os.path.join(EDGE_DIR, "*_edge_list.csv")))
    if not edge_files:
        print(f"No edge lists found in: {EDGE_DIR}")
        return

    print(f"Found {len(edge_files)} edge lists.")
    print(f"PNG output dir: {PNG_DIR}")

    successes = 0
    failures = 0

    for edge_csv in edge_files:
        prefix = prefix_from_edge_path(edge_csv)
        node_csv = os.path.join(NODE_DIR, f"{prefix}_node_attributes.csv")
        out_png = os.path.join(PNG_DIR, f"{prefix}.png")

        print(f"\n--- {prefix} ---")
        print(f"Edge: {edge_csv}")
        print(f"Node: {node_csv}")
        print(f"PNG:  {out_png}")

        if not os.path.exists(node_csv):
            print(f"SKIP: missing node attributes file: {node_csv}")
            failures += 1
            continue

        try:
            net_suid = import_network_from_csv(edge_csv, title=prefix)
            import_node_table(node_csv, net_suid)

            apply_style_to_network(net_suid)
            layout_network(net_suid)
            time.sleep(0.3)

            export_png(net_suid, out_png)
            print("OK: exported PNG")
            successes += 1

            delete_network(net_suid)

        except Exception as e:
            failures += 1
            print("ERROR:", e)
            traceback.print_exc()
            continue

    print("\nDone.")
    print(f"Successes: {successes}")
    print(f"Failures:  {failures}")

if __name__ == "__main__":
    main()
