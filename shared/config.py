from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = Path(os.getenv("DEAPP_DATA_ROOT", PROJECT_ROOT / "data"))
REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
RUN_INLINE = os.getenv("RUN_INLINE", "false").lower() in {"1", "true", "yes"}

# Optional S3-compatible object storage settings.
# When S3_BUCKET is set, uploads/artifacts/metadata are stored in object storage.
S3_BUCKET = os.getenv("S3_BUCKET", "").strip()
S3_REGION = os.getenv("S3_REGION", "us-east-1")
S3_ENDPOINT_URL = os.getenv("S3_ENDPOINT_URL", "").strip() or None
S3_ACCESS_KEY_ID = os.getenv("S3_ACCESS_KEY_ID", "").strip() or None
S3_SECRET_ACCESS_KEY = os.getenv("S3_SECRET_ACCESS_KEY", "").strip() or None
S3_PREFIX = os.getenv("S3_PREFIX", "deapp").strip().strip("/")
S3_PRESIGN_EXPIRES = int(os.getenv("S3_PRESIGN_EXPIRES", "3600"))
