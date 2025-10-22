import os
import json
import pandas as pd

def load_config(config_path, path_vars):
    """Improved config loader with path replacement and fallback resolution"""
    try:
        with open(config_path) as f:
            config = json.load(f)

        if not config:
            raise ValueError("Config file is empty")

        # Recursively replace placeholders in the config
        def replace_placeholders(obj):
            if isinstance(obj, dict):
                return {k: replace_placeholders(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [replace_placeholders(item) for item in obj]
            elif isinstance(obj, str):
                try:
                    return obj.format(**path_vars)
                except KeyError:
                    return obj  # Return unchanged if placeholder doesn't exist
            return obj

        config = replace_placeholders(config)

        # Handle BASE_SS fallback logic if the path doesn't exist
        if "paths" in config and "BASE_SS" in config["paths"]:
            base_ss_path = config["paths"]["BASE_SS"]
            if not os.path.exists(base_ss_path):
                # Check fallbacks if they exist
                fallbacks = config["paths"].get("_BASE_SS_fallbacks", [])
                for fallback in fallbacks:
                    print(fallback)
                    print(os.path.exists(fallback))
                    if os.path.exists(fallback):
                        config["paths"]["BASE_SS"] = fallback
                        break
                else:
                    print(f"Warning: None of the BASE_SS paths exist ({base_ss_path} or fallbacks)")

        return config

    except FileNotFoundError:
        print(f"Error: Config file not found at {config_path}")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON in config file at {config_path}")
        return None
    except Exception as e:
        print(f"Error loading config: {str(e)}")
        return None

def extract_atlas_definitions(config):
    """Extract all atlas definitions as DataFrames"""
    atlas_dfs = {}
    for atlas_name in config['atlas_definitions']:
        atlas_def = config['atlas_definitions'][atlas_name]
        regions = atlas_def['regions']

        if atlas_def['type'] == 'bilateral':
            df = pd.DataFrame([{
                'region': r['region'],
                'label': r['label'],
                'hemisphere': 'both'
            } for r in regions])
        else:
            df = pd.DataFrame([{
                'region': r['region'],
                'label': r['label'],
                'hemisphere': 'left' if r['region'].startswith('L_') else 'right'
            } for r in regions])

        atlas_dfs[atlas_name] = df

    return atlas_dfs