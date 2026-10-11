"""Structured validation reports. Validators never fix input silently; they say what and where."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Issue:
    """One problem found in user input.

    ``where`` locates it: an entity name and field (``"spline-1: points[3]"``) or a file and
    line (``"naca0012.dat:14"``).
    """

    where: str
    message: str

    def __str__(self) -> str:
        return f"{self.where}: {self.message}"


class ValidationError(ValueError):
    """Raised with one or more issues; ``str()`` lists them all."""

    def __init__(self, issues: list[Issue]):
        self.issues = list(issues)
        super().__init__("\n".join(str(issue) for issue in self.issues))
