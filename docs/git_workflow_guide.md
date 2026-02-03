# Git Workflow Guide for Lysis Project

## Philosophy

Keep it simple and pragmatic:
- **Branches** isolate work and make it safe to experiment
- **Regular commits** create save points you can return to
- **Descriptive messages** help you remember what you did
- **Don't stress** about perfect history - progress > perfection

## Branch Naming Convention

Use descriptive names with a category prefix:

```bash
feature/task-name       # New features or functionality
fix/bug-description     # Bug fixes
docs/what-documenting   # Documentation updates
refactor/what-refactoring  # Code restructuring
experiment/what-trying  # Experimental changes
```

**Examples:**
```bash
feature/datastore-implementation
fix/hdf5-attribute-writing
docs/data-handling-completion-plan
refactor/extract-constants
```

## When to Create a Branch

Create a new branch when you're about to:

1. **Start a new feature** (Task #1, #5, etc.)
   ```bash
   git checkout -b feature/datastore-core-methods
   ```

2. **Fix a bug** (Task #3, #4)
   ```bash
   git checkout -b fix/hdf5-attribute-slice
   ```

3. **Try something experimental**
   ```bash
   git checkout -b experiment/alternative-converter-approach
   ```

4. **Work on something that might take multiple sessions**
   ```bash
   git checkout -b feature/macroscale-converters
   ```

**Stay on main for:**
- Quick typo fixes
- Updating documentation (unless extensive)
- Small formatting changes

## When to Commit

Commit after completing a **logical unit of work**:

### ✅ Good times to commit:

- ✅ Added a new function/method
- ✅ Fixed a bug
- ✅ Updated documentation
- ✅ Refactored a section
- ✅ Completed a TODO item
- ✅ Tests pass
- ✅ End of work session (even if incomplete)

### ❌ Don't wait for:

- ❌ Everything to be perfect
- ❌ All tests to be written
- ❌ Documentation to be complete
- ❌ End of entire feature

### Commit Message Template:

```bash
# Short (50 chars or less) summary in imperative mood
Add _not_implemented helper function to dataconvert

# Blank line, then details if needed
Creates converter stubs that raise NotImplementedError with clear
messages showing dataset name and spec versions. Added stubs for
20+ missing macroscale converters.

# Optional: Reference issues/tasks
Relates to Task #15, #16
```

**Quick commits** (one-liners are fine):
```bash
git commit -m "Fix HDF5 attribute writing bug"
git commit -m "Remove debug print statement"
git commit -m "Add docstrings to parameters.py"
```

## Commit Frequency Guidelines

**Aim for commits every 20-40 minutes of focused work**, or:

- After ~50-100 lines of new code
- After each method/function completed
- Before switching tasks
- Before trying something risky
- Before taking a break

**Rule of thumb:** If you'd be upset to lose the work you just did, commit it!

## When to Merge

Merge a feature branch when:

1. **Feature is complete** (passes your own testing)
2. **Code works** (doesn't break existing functionality)
3. **Ready to move on** (don't leave stale branches)

### Simple Merge Workflow:

```bash
# 1. Make sure your branch is up to date and committed
git status
git add .
git commit -m "Complete feature X"

# 2. Switch to main
git checkout main

# 3. Pull any changes (if working with others)
git pull

# 4. Merge your feature branch
git merge feature/your-feature-name

# 5. Push to remote
git push

# 6. Delete the feature branch (cleanup)
git branch -d feature/your-feature-name
```

### When NOT to merge yet:

- ❌ Feature is half-done and you need a break → commit and push the branch
- ❌ Tests are failing → fix first, then merge
- ❌ You want to try an alternative approach → keep branch, create another

## Example Workflow for Task #1 (DataStore Implementation)

```bash
# 1. Create branch for the feature
git checkout -b feature/datastore-core-methods

# 2. Work on __getattr__ method
#    ... write code ...
git add src/python/lysis/util/datastore.py
git commit -m "Implement DataStore.__getattr__ method"

# 3. Work on __setattr__ method
#    ... write code ...
git add src/python/lysis/util/datastore.py
git commit -m "Implement DataStore.__setattr__ method"

# 4. Work on status() method
#    ... write code ...
git add src/python/lysis/util/datastore.py
git commit -m "Implement DataStore.status method"

# 5. Add tests
#    ... write tests ...
git add tests/test_datastore.py
git commit -m "Add tests for DataStore core methods"

# 6. Feature complete! Merge to main
git checkout main
git merge feature/datastore-core-methods
git push

# 7. Cleanup
git branch -d feature/datastore-core-methods
```

## Handling Mistakes

### Forgot to create a branch?

```bash
# Create branch from current state
git checkout -b feature/forgot-to-branch
# Your changes are now on the new branch!
```

### Want to undo last commit?

```bash
# Keep changes, undo commit
git reset --soft HEAD~1

# Discard changes and commit
git reset --hard HEAD~1  # ⚠️ Destructive!
```

### Made commits on main instead of feature branch?

```bash
# Create branch with current state
git checkout -b feature/should-have-branched

# Reset main to before your commits
git checkout main
git reset --hard origin/main

# Your work is safe on the feature branch
git checkout feature/should-have-branched
```

## Quick Reference

### Most Common Commands

```bash
# Check status
git status

# Create and switch to branch
git checkout -b feature/branch-name

# Stage changes
git add <file>              # Specific file
git add .                   # All changes

# Commit
git commit -m "Message"

# Switch branches
git checkout main
git checkout feature/branch-name

# Merge (from main)
git merge feature/branch-name

# Push
git push
git push -u origin feature/branch-name  # First time pushing branch

# List branches
git branch                  # Local branches
git branch -a              # All branches

# Delete branch
git branch -d feature/branch-name
```

## Current Work Status

Use this to track what you're working on:

- **Current branch**: _____________
- **Current task**: _____________
- **Last commit**: _____________
- **Ready to merge?**: Yes / No / Not sure

---

## Reminders from Claude

I'll remind you to:

- 🔀 **Branch** when starting a new task
- 💾 **Commit** after completing logical units of work
- 🔄 **Merge** when a feature is complete
- 🧹 **Clean up** old branches after merging

Don't hesitate to ask me:
- "Should I branch for this?"
- "Is this a good time to commit?"
- "Should I merge this now?"
- "What should my commit message be?"
