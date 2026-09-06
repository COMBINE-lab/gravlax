"""Bound report size without weakening the complete-result equality checks."""
import copy
import hashlib
import json


def compact_summary(summary):
    result = copy.deepcopy(summary)
    if result and len(result.get('population_by_group', [])) > 100:
        groups = result.pop('population_by_group')
        result['population_by_group_digest'] = {
            'count': len(groups),
            'sha256': hashlib.sha256(json.dumps(groups, sort_keys=True).encode()).hexdigest(),
            'note': 'Complete populations compared before reporting; large group list stored as digest.',
        }
    return result
