"""Cross-cutting services shared by strategies: RNG handling and recording.
These are not strategies (no move/bind/unbind polymorphism) -- they're
plain dependencies that get constructed once in factory.py and passed
into whichever strategies need them.
"""
