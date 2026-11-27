"""FLAAP launcher entry point (AI-CAE4ALL style).

    python FLAAP_main.py --config configs/flutter/ex1/config_goland_pk.txt
    python FLAAP_main.py --list-models
    python FLAAP_main.py --describe flutter
"""
import sys

from flaap.launcher.main import main

if __name__ == "__main__":
    sys.exit(main())
