"""Root conftest: add skill directory to sys.path so tests find scripts/."""
import sys
from pathlib import Path

skill_dir = Path(__file__).parent / "skills" / "hirn-literature-retrieve"
sys.path.insert(0, str(skill_dir))
