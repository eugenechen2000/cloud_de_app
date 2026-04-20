from __future__ import annotations

import redis
from rq import Queue

from shared.config import REDIS_URL


def get_queue() -> Queue:
    conn = redis.from_url(REDIS_URL)
    return Queue("de-jobs", connection=conn)
