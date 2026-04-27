"""
Schema validation utilities for industrial lignin JSON input.
"""

from __future__ import annotations

import json
import os
from typing import Any

import jsonschema


class InputSchemaValidator:
    """
    Validate input dictionaries/files against `lignin_info_schema.json`.
    """

    def __init__(self, schema_path: str | None = None):
        if schema_path is None:
            schema_path = os.path.join(os.path.dirname(__file__), "lignin_info_schema.json")
        self.schema_path = schema_path
        with open(self.schema_path, "r") as f:
            self.schema = json.load(f)

    def validate_dict(self, data: dict[str, Any]) -> None:
        jsonschema.validate(instance=data, schema=self.schema)

    def validate_file(self, json_path: str) -> dict[str, Any]:
        with open(json_path, "r") as f:
            data = json.load(f)
        self.validate_dict(data)
        return data
