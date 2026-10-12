"""The topology-editing command: split / join / O-grid operations (opengrid.model.topology_edit)
as one undo step, journaled as the high-level call (the operations are deterministic, so replay
creates the same entities with the same names)."""

from opengrid.commands.base import Command
from opengrid.commands.geometry import format_value
from opengrid.model.store import Store


class TopologyEdit(Command):
    def __init__(self, method: str, function, args: tuple, kwargs: dict, label: str):
        self.method = method
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.label = label
        self.warnings: list[str] = []
        self._script = ""

    def execute(self, store: Store) -> list[str]:
        names, warnings = self.function(store, *self.args, **self.kwargs)
        self.warnings = list(warnings)
        parts = [format_value(a) for a in self.args]
        parts += [f"{k}={format_value(v)}" for k, v in self.kwargs.items()]
        self._script = f"og.{self.method}({', '.join(parts)})"
        return names

    def script(self) -> str:
        return self._script
