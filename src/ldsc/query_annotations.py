"""Internal records for batch query-annotation status and provenance."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any


@dataclass(frozen=True)
class QueryAnnotationStatus:
    """Status of one requested BED or gene-list query annotation."""

    query: str
    source: str
    input_type: str
    status: str
    reason: str = ""
    n_annotation_snps: float | None = None
    details: str | None = None

    def updated(self, **changes: Any) -> "QueryAnnotationStatus":
        """Return a copy with final workflow status fields replaced."""
        return replace(self, **changes)

    def as_dict(self) -> dict[str, Any]:
        """Return the fixed query-status manifest row schema."""
        return {
            "query": self.query,
            "source": self.source,
            "input_type": self.input_type,
            "status": self.status,
            "reason": self.reason,
            "n_annotation_snps": self.n_annotation_snps,
            "details": self.details,
        }
