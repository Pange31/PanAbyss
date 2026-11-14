import sqlite3
import os
import uuid
import logging

from datetime import datetime, timezone
from config import *

logger = logging.getLogger("panabyss_logger")

DB_PATH = os.path.abspath("./sqlite")
DB_FILENAME = os.path.join(DB_PATH, "management_jobs.db")


def now_utc():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def get_management_connection():
    os.makedirs(DB_PATH, exist_ok=True)

    return sqlite3.connect(
        DB_FILENAME,
        timeout=30,
        check_same_thread=False
    )


def init_management_db():
    conn = get_management_connection()

    try:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS management_jobs (
                job_id TEXT PRIMARY KEY,
                status TEXT NOT NULL,
                created_at DATETIME,
                updated_at DATETIME
            )
        """)

        conn.commit()
        logger.debug(f"Database initialized at {DB_FILENAME}")

    finally:
        conn.close()


def create_index_job():
    job_id = str(uuid.uuid4())
    conn = get_management_connection()

    try:
        timestamp = now_utc()

        conn.execute("""
            INSERT INTO management_jobs
            (job_id, status, created_at, updated_at)
            VALUES (?, 'running', ?, ?)
        """, (job_id, timestamp, timestamp))

        conn.commit()

        logger.info(f"Index creation job created: {job_id}")

        return job_id

    except Exception:
        conn.rollback()
        logger.exception("Unable to create index job")
        raise

    finally:
        conn.close()


def request_index_cancellation(job_id):
    conn = get_management_connection()

    try:
        cursor = conn.execute("""
            UPDATE management_jobs
            SET status = 'cancel_requested',
                updated_at = ?
            WHERE job_id = ?
            AND status = 'running'
        """, (now_utc(), job_id))

        conn.commit()

        logger.info(
            f"Cancellation requested for job {job_id}. "
            f"Updated rows: {cursor.rowcount}"
        )

        return cursor.rowcount > 0

    except Exception:
        conn.rollback()
        logger.exception("Unable to request index cancellation")
        raise

    finally:
        conn.close()


def is_index_cancellation_requested(job_id):
    conn = get_management_connection()

    try:
        result = conn.execute("""
            SELECT status
            FROM management_jobs
            WHERE job_id = ?
        """, (job_id,)).fetchone()

        return result is not None and result[0] == "cancel_requested"

    finally:
        conn.close()


def update_index_job_status(job_id, status):
    conn = get_management_connection()

    try:
        cursor = conn.execute("""
            UPDATE management_jobs
            SET status = ?, updated_at = ?
            WHERE job_id = ?
        """, (status, now_utc(), job_id))

        conn.commit()

        logger.info(
            f"Job {job_id} status updated to {status}. "
            f"Updated rows: {cursor.rowcount}"
        )

    except Exception:
        conn.rollback()
        logger.exception("Unable to update index job status")
        raise

    finally:
        conn.close()

def has_active_index_job():
    conn = get_management_connection()

    try:
        result = conn.execute("""
            SELECT COUNT(*)
            FROM management_jobs
            WHERE status IN ('running', 'cancel_requested')
        """).fetchone()

        return result[0] > 0

    finally:
        conn.close()

#This function will delete all the table
def drop_management_table():
    conn = get_management_connection()
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS management_jobs")

    conn.commit()
    conn.close()