#This file allows to defines functions used for retrocompatibility


import os
import platform
import shutil
import logging
from pathlib import Path
import re


logger = logging.getLogger("panabyss_logger")

def check_updates():
    #check required update scripts to update from <= 1.5.0 to > 1.6.0 versions
    migrate_neo4j_directory_structure()

# This file defines functions used for retrocompatibility



logger = logging.getLogger("panabyss_logger")

"""
Check and apply required updates for retrocompatibility.
"""
def check_updates():
    migrate_neo4j_directory_structure()
    deactivate_neo4j_auth()

"""
    For migration from <= 1.5.0 to >= 1.6.0 versions
    Migrate the old Neo4j directory structure.
    
    Old structure:
        data/
            data/
            logs/
            plugins/
    
    New structure:
        data/
            database/
                data/
                logs/
                plugins/
    
    The migration supports partially migrated installations.
    
    If at least one directory is migrated:
      - the old Docker Compose file is deleted;
      - all files and directories under data/ are chmod 777.
    
    Compatible with Linux, macOS and Windows.
"""
def migrate_neo4j_directory_structure():
    # ------------------------------------------------------------------
    # Determine application directories
    # ------------------------------------------------------------------

    app_dir = Path(__file__).resolve().parents[1]
    base_dir = app_dir / "data"
    database_dir = base_dir / "database"

    docker_compose_file = (
            app_dir
            / "docker_management"
            / "docker-compose.yml"
    )

    # ------------------------------------------------------------------
    # OLD DIRECTORIES
    # ------------------------------------------------------------------

    old_data_dir = base_dir / "data"
    old_logs_dir = base_dir / "logs"
    old_plugins_dir = base_dir / "plugins"

    if not (
            old_data_dir.is_dir()
            or old_logs_dir.is_dir()
            or old_plugins_dir.is_dir()
    ):
        logger.debug(
            "✅ No old Neo4j directories found. "
            "No migration required."
        )
        return

    logger.info(
        "🔄 Old Neo4j directory structure detected. "
        "Starting migration..."
    )

    # ------------------------------------------------------------------
    # Define migrations
    # ------------------------------------------------------------------

    migrations = [
        (
            old_data_dir,
            database_dir / "data",
        ),
        (
            old_logs_dir,
            database_dir / "logs",
        ),
        (
            old_plugins_dir,
            database_dir / "plugins",
        ),
    ]

    # ------------------------------------------------------------------
    # Create database directory only when migration is actually needed
    # ------------------------------------------------------------------

    try:
        database_dir.mkdir(
            parents=True,
            exist_ok=True,
        )
    except OSError as exc:
        logger.error(
            f"❌ Cannot create database directory "
            f"{database_dir}: {exc}"
        )
        raise

    migration_done = False

    # ------------------------------------------------------------------
    # Process each old directory independently
    # ------------------------------------------------------------------

    for source, destination in migrations:

        # Source does not exist.
        # This is normal for a partially migrated installation.
        if not source.is_dir():
            logger.debug(
                f"⏭️ Source does not exist: {source}"
            )
            continue

        # --------------------------------------------------------------
        # Destination already exists.
        #
        # We keep the destination and delete the old source.
        # --------------------------------------------------------------

        if destination.exists():

            logger.info(
                f"⚠️ Destination already exists: {destination}"
            )

            logger.info(
                f"🗑️ Removing old directory: {source}"
            )

            try:
                shutil.rmtree(source)
            except OSError as exc:
                logger.error(
                    f"❌ Cannot remove old directory "
                    f"{source}: {exc}"
                )
                raise

            logger.info(
                f"✅ Old directory removed: {source}"
            )

            migration_done = True
            continue

        # --------------------------------------------------------------
        # Destination does not exist.
        #
        # shutil.move() does not normally copy the content when source
        # and destination are on the same filesystem.
        # --------------------------------------------------------------

        logger.info(
            f"🔄 Moving:\n"
            f"   {source}\n"
            f"   → {destination}"
        )

        try:
            destination.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            shutil.move(
                str(source),
                str(destination),
            )

        except OSError as exc:
            logger.error(
                f"❌ Cannot move:\n"
                f"   {source}\n"
                f"   → {destination}\n"
                f"   {exc}"
            )
            raise

        logger.info(
            f"✅ Migration completed: {destination}"
        )

        migration_done = True

    # ------------------------------------------------------------------
    # If something changed, delete Docker Compose file.
    # ------------------------------------------------------------------

    if not migration_done:
        return

    if docker_compose_file.exists():

        logger.info(
            f"🗑️ Deleting Docker Compose file: "
            f"{docker_compose_file}"
        )

        try:
            docker_compose_file.unlink()
        except OSError as exc:
            logger.error(
                f"❌ Cannot delete Docker Compose file: "
                f"{docker_compose_file}: {exc}"
            )
            raise

    # ------------------------------------------------------------------
    # Set permissions to 777 recursively.
    #
    # This is intentionally done ONLY after a migration.
    # Therefore a normal startup on an already migrated database
    # does not scan the entire database.
    # ------------------------------------------------------------------

    logger.info(
        f"🔐 Setting permissions to 777 recursively: "
        f"{base_dir}"
    )

    for root, dirs, files in os.walk(base_dir):

        root_path = Path(root)

        try:
            os.chmod(root_path, 0o777)
        except OSError as exc:
            logger.warning(
                f"⚠️ Cannot chmod {root_path}: {exc}"
            )

        for name in dirs:

            path = root_path / name

            try:
                os.chmod(path, 0o777)
            except OSError as exc:
                logger.warning(
                    f"⚠️ Cannot chmod {path}: {exc}"
                )

        for name in files:

            path = root_path / name

            try:
                os.chmod(path, 0o777)
            except OSError as exc:
                logger.warning(
                    f"⚠️ Cannot chmod {path}: {exc}"
                )

    logger.info(
        "✅ Neo4j directory migration completed."
    )

"""
    For migration from <= 1.5.0 to >= 1.6.0 versions
    deactivate neo4j authentification
    dbms.security.auth_enabled is set to false
"""
def deactivate_neo4j_auth():
    app_dir = Path(__file__).resolve().parents[1]
    neo4j_conf_file = app_dir / "data" / "conf" / "neo4j.conf"
    conf_file = Path(neo4j_conf_file)
    parameter = "dbms.security.auth_enabled"

    if not conf_file.exists():
        logger.error(f"❌ Neo4j config file not found: {conf_file}")
        return False

    try:
        content = conf_file.read_text(encoding="utf-8")
        pattern = re.compile(
            rf"^(\s*#?\s*{re.escape(parameter)}\s*=).*$",
            re.MULTILINE
        )

        if pattern.search(content):
            content = pattern.sub(
                rf"\1false",
                content
            )
        else:
            if content and not content.endswith("\n"):
                content += "\n"

            content += f"{parameter}=false\n"

        conf_file.write_text(content, encoding="utf-8")

        logger.info(
            f"🔓 Neo4j authentication disabled in {conf_file}"
        )

        return True

    except Exception as e:
        logger.error(
            f"❌ Failed to disable Neo4j authentication: {e}"
        )
        return False
