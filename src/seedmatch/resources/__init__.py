"""The registry of the upstream resources that the `seedmatch-fetch` command downloads."""

from seedmatch.resources.registry import (
    REGISTRY,
    STRING_VERSION,
    Entry,
    Group,
    Resource,
    lookup,
)

__all__ = ["REGISTRY", "STRING_VERSION", "Entry", "Group", "Resource", "lookup"]
