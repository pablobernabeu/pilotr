<!-- Thank you for contributing to pilotr. -->

**What does this change do?**
A short description, and a link to the issue it addresses if there is one.

**Checklist**

R package (`r/pilotr/`):

- [ ] `devtools::document()` has been run if the roxygen comments changed.
- [ ] `devtools::test()` passes, and new behaviour has tests.
- [ ] `devtools::check()` is clean.

Python package (`python/`):

- [ ] `pytest` passes, and new behaviour has tests.
- [ ] The docs build (`mkdocs build --strict`).

Both:

- [ ] Changes to the generative core land in both implementations, and the
      `tools/parity/` harness passes, as the contributing guide sets out.
- [ ] Documentation and `NEWS.md` / `CHANGELOG.md` are updated where relevant.
- [ ] No secret appears anywhere in the diff.
