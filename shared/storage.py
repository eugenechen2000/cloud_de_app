from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import boto3
from botocore.client import BaseClient

from shared.config import (
    DATA_ROOT,
    S3_ACCESS_KEY_ID,
    S3_BUCKET,
    S3_ENDPOINT_URL,
    S3_PREFIX,
    S3_PRESIGN_EXPIRES,
    S3_REGION,
    S3_SECRET_ACCESS_KEY,
)


@dataclass
class StorageBackend:
    mode: str  # "local" or "s3"
    bucket: str | None
    prefix: str
    root: Path
    s3: BaseClient | None

    def key(self, rel_path: str) -> str:
        rel = rel_path.lstrip("/")
        return f"{self.prefix}/{rel}" if self.prefix else rel

    def put_bytes(self, rel_path: str, data: bytes, content_type: str | None = None) -> None:
        if self.mode == "local":
            path = self.root / rel_path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            return
        extra = {}
        if content_type:
            extra["ContentType"] = content_type
        self.s3.put_object(Bucket=self.bucket, Key=self.key(rel_path), Body=data, **extra)

    def get_bytes(self, rel_path: str) -> bytes:
        if self.mode == "local":
            return (self.root / rel_path).read_bytes()
        obj = self.s3.get_object(Bucket=self.bucket, Key=self.key(rel_path))
        return obj["Body"].read()

    def exists(self, rel_path: str) -> bool:
        if self.mode == "local":
            return (self.root / rel_path).exists()
        try:
            self.s3.head_object(Bucket=self.bucket, Key=self.key(rel_path))
            return True
        except Exception:
            return False

    def list_files(self, rel_prefix: str) -> list[str]:
        rel_prefix = rel_prefix.rstrip("/") + "/"
        if self.mode == "local":
            base = self.root / rel_prefix
            if not base.exists():
                return []
            return sorted([p.name for p in base.glob("*") if p.is_file()])

        out: list[str] = []
        paginator = self.s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=self.bucket, Prefix=self.key(rel_prefix)):
            for obj in page.get("Contents", []):
                k = obj["Key"]
                suffix = k.split("/", 1)[1] if "/" in k else k
                if not suffix.startswith(rel_prefix):
                    continue
                name = suffix[len(rel_prefix):]
                if name and "/" not in name:
                    out.append(name)
        return sorted(set(out))

    def file_size(self, rel_path: str) -> int:
        if self.mode == "local":
            return (self.root / rel_path).stat().st_size
        head = self.s3.head_object(Bucket=self.bucket, Key=self.key(rel_path))
        return int(head["ContentLength"])

    def download_to(self, rel_path: str, local_path: Path) -> None:
        local_path.parent.mkdir(parents=True, exist_ok=True)
        if self.mode == "local":
            local_path.write_bytes((self.root / rel_path).read_bytes())
            return
        self.s3.download_file(self.bucket, self.key(rel_path), str(local_path))

    def upload_from(self, local_path: Path, rel_path: str, content_type: str | None = None) -> None:
        if self.mode == "local":
            dst = self.root / rel_path
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(local_path.read_bytes())
            return
        extra = {}
        if content_type:
            extra["ExtraArgs"] = {"ContentType": content_type}
        self.s3.upload_file(str(local_path), self.bucket, self.key(rel_path), **extra)

    def presigned_url(self, rel_path: str, expires_seconds: int | None = None) -> str | None:
        if self.mode == "local":
            return None
        exp = expires_seconds or S3_PRESIGN_EXPIRES
        return self.s3.generate_presigned_url(
            ClientMethod="get_object",
            Params={"Bucket": self.bucket, "Key": self.key(rel_path)},
            ExpiresIn=exp,
        )


def get_storage() -> StorageBackend:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    if not S3_BUCKET:
        return StorageBackend(mode="local", bucket=None, prefix="", root=DATA_ROOT, s3=None)

    session = boto3.session.Session(
        aws_access_key_id=S3_ACCESS_KEY_ID,
        aws_secret_access_key=S3_SECRET_ACCESS_KEY,
        region_name=S3_REGION,
    )
    s3 = session.client("s3", endpoint_url=S3_ENDPOINT_URL)
    return StorageBackend(mode="s3", bucket=S3_BUCKET, prefix=S3_PREFIX, root=DATA_ROOT, s3=s3)
