from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import logging
import threading
import logging
from config import get_conf
import time

logger = logging.getLogger("panabyss_logger")
#main driver
_DRIVER = None
#To avoid collision between thread when getting the driver
_LOCK = threading.Lock()

CONF = None

def get_driver(max_retries=10, retry_delay=3):
    global _DRIVER
    global CONF

    # 🔥 fast path (sans lock)
    if _DRIVER is not None:
        return _DRIVER

    if not CONF:
        CONF = get_conf()

    if "container_name" not in CONF or not CONF["container_name"]:
        return None

    DB_URL = CONF["DB_URL"]
    AUTH = CONF["AUTH"]

    # 🔥 PROTECTION MULTI-THREAD
    with _LOCK:

        # re-check après lock (IMPORTANT)
        if _DRIVER is not None:
            return _DRIVER

        for attempt in range(1, max_retries + 1):
            try:
                driver = GraphDatabase.driver(
                    DB_URL,
                    auth=AUTH,
                    max_connection_lifetime=3600,
                    liveness_check_timeout=10,
                    max_connection_pool_size=100,
                )

                with driver.session() as session:
                    session.run("RETURN 1").consume()

                _DRIVER = driver
                logger.info("✅ Neo4j driver initialized and ready")
                return _DRIVER

            except Exception as e:
                logger.warning(
                    f"⏳ Neo4j not ready (attempt {attempt}/{max_retries}): {e}"
                )

                try:
                    driver.close()
                except:
                    pass

                time.sleep(retry_delay)

    logger.error("❌ Failed to connect to Neo4j after retries")
    return None

def reset_driver():
    global _DRIVER
    with _LOCK:
        if _DRIVER:
            try:
                _DRIVER.close()
            except:
                pass
        _DRIVER = None

#Create a  Neo4j driver instance (NOT global),
#intended for heavy jobs / isolated workloads.
def get_scoped_driver(max_retries=10, retry_delay=3):
    if not CONF:
        conf = get_conf()
    else:
        conf = CONF

    if "container_name" not in conf or not conf["container_name"]:
        return None

    DB_URL = conf["DB_URL"]
    AUTH = conf["AUTH"]

    for attempt in range(1, max_retries + 1):
        driver = None
        try:
            driver = GraphDatabase.driver(
                DB_URL,
                auth=AUTH,
                max_connection_lifetime=3600,
                liveness_check_timeout=10,
                max_connection_pool_size=20,   # 🔥 plus petit conseillé par job
                connection_acquisition_timeout=30,
            )

            # 🔥 test immédiat
            with driver.session() as session:
                session.run("RETURN 1").consume()

            logger.info(f"✅ Scoped Neo4j driver ready (attempt {attempt})")
            return driver

        except Exception as e:
            logger.warning(
                f"⏳ Scoped driver init failed (attempt {attempt}/{max_retries}): {e}"
            )

            try:
                if driver:
                    driver.close()
            except Exception:
                pass

            time.sleep(retry_delay)

    logger.error("❌ Failed to create scoped Neo4j driver")
    return None