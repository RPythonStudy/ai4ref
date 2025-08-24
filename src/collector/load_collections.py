
import json
from pathlib import Path
from common.database import get_db_connection

def load_search_collections():
    json_path = Path("config/search_collections.json")
    with open(json_path, "r", encoding="utf-8") as f:
        collections = json.load(f)["search_collections"]
    conn = get_db_connection()
    cur = conn.cursor()
    for table in ["paper_collection", "collection_pmid", "unique_pmid", "filtered_pmid", "collections"]:
        cur.execute(f"DELETE FROM {table}")
    for c in collections:
        cur.execute(
            """
            INSERT INTO collections (name, search_term, retmax, enabled, description, created_at, updated_at)
            VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
            """,
            (
                c.get("collection"),
                c.get("term"),
                c.get("retmax", 10000),
                c.get("enabled", True),
                c.get("description", "")
            )
        )
    conn.commit()
    cur.close()
    conn.close()

if __name__ == "__main__":
    load_search_collections()
