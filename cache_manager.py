from diskcache import Cache
import uuid
from datetime import datetime
import logging

logger = logging.getLogger("panabyss_logger")

#create a cache limited to 500 Mo
nodes_cache = Cache("./cache", size_limit=500 * 1024 * 1024)

def get_session_cache(nodes_cache_id):
    cached = nodes_cache.get(nodes_cache_id)
    logger.debug(f"cached id : {nodes_cache_id}")
    if cached is None:
        cached = {
            "nodes": {},
            "zoom": {},
            "min_node_size": None,
            "last_update": datetime.now().isoformat()
        }

        nodes_cache.set(
            nodes_cache_id,
            cached,
            expire=8 * 3600
        )

    return cached