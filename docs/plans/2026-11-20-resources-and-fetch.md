# Resources and fetch: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the registry of the upstream files, the download module, and the command `seedmatch-fetch`. A machine with network access uses them to prepare each file that the sandbox reads from a path.

**Architecture:** `seedmatch.resources` holds the registry as frozen data, with no I/O. `seedmatch.fetch` streams a registry file to `<dest>/<name>.part`, makes sure that the size and the SHA-256 digest agree, and then renames the file. The HTTP client is a function in `FetchOpts`, thus the tests serve fixture bytes from memory. `seedmatch.fetch.cli` gives the command.

**Tech Stack:** Python 3.12 standard library only: `urllib.request`, `hashlib`, `argparse`, `os`, `pathlib`.

**Spec:** `docs/design.md`, sections "Architecture", "Resources and fetch", and "Tests" ("Core without network", "Fetch").

**Later plans (not in this plan):** the readers of miRDB, TRRUST, and STRING, and the BioMart map of the STRING protein IDs. These come with the analyses that use them.

## Global Constraints

- No new runtime dependency. The fetch module uses `urllib.request` from the standard library.
- Only `seedmatch.fetch` imports a network module. The ruff rule `TID251` bans `urllib.request`, `http.client`, `socket`, `ssl`, `urllib3`, `requests`, `httpx`, and `seedmatch.fetch` in each other module. `import seedmatch` does not import `seedmatch.fetch`.
- The tests never use the network (`pytest-socket`, `--disable-socket`). The tests give `FetchOpts` an HTTP client that serves fixture bytes from memory.
- Errors derive from `seedmatch.errors.SeedMatchError`. Invalid input raises `InputError`. A failed download raises `FetchError`, a new subclass of `SeedMatchError`. Each message names the bad value and the expected value.
- A function writes a file only into its explicit `dest` directory.
- Documents use ASD-STE100. Code comments obey the admission test of the project rules. The word "seam" does not occur.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits.

## The registry

The controller streamed each file on 2026-11-20 and measured the size and the SHA-256 digest (`tmp/claude/resources/hash_remote.py`). No URL redirected. The file name of each resource is the last part of its URL path.

Reference builds (group `reference`). The URL is `location_prefix + rdatapath` of `https://annotationhub.bioconductor.org/ahid/<id>`, with the prefix `https://mghp.osn.xsede.org/bir190004-bucket01/AnnotationHub/`. The table `ah.df` of `R/load_annotations.R:152-154` gives the IDs.

| Key | AnnotationHub ID | rdatapath | Size | SHA-256 |
|---|---|---|---|---|
| `hg38` | AH106283 | `ensembl/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.2bit` | 774952244 | `6464ece374014a850de7da50bfaab5083b73197bdeb1393731065b6ab7560be6` |
| `hg38` | AH109606 | `AHEnsDbs/v109/EnsDb.Hsapiens.v109.sqlite` | 445898752 | `a3d9a352b0bc07ac6f27e89df31308ff8ead36380608644cc5939de3212b8253` |
| `mm39` | AH106440 | `ensembl/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.2bit` | 682060023 | `eac802b955d7885f0eb620ba5d49f119397bf5c630593e20d16f0f5bf3d4e8e9` |
| `mm39` | AH109655 | `AHEnsDbs/v109/EnsDb.Mmusculus.v109.sqlite` | 269914112 | `6d2efcbfd5d3cb14ba2037d67583e1d4027b1eba81d489a9a4938f948107a7c2` |
| `mm10` | AH88475 | `ensembl/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.2bit` | 682724912 | `b12e356f498d74a8af43e374d3240ab058d1d6b9f76baefb277d5edb469f488a` |
| `mm10` | AH89211 | `AHEnsDbs/v102/EnsDb.Mmusculus.v102.sqlite` | 265134080 | `b8b915e96131fadf494d3f7a79e9e7b31a54aa002cc60e4ba2278797a2ce0c08` |
| `rnor6` | AH93578 | `ensembl/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.2bit` | 718742859 | `e5eb2cd9a17d8f52c33d07bb423085addbc74a3c5ec3b1b2a8f2c6eb0ef6a4ee` |
| `rnor6` | AH95846 | `AHEnsDbs/v104/EnsDb.Rnorvegicus.v104.sqlite` | 117067776 | `174312139350c7cb86af5fd82c9758c98c4fe43d6f75ca21797f2c154fa1bdae` |
| `rnor7` | AH106786 | `ensembl/release-107/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit` | 661989722 | `73eacccd66a8843cced837d0cb0a272b42b7214b4a7c62a3c48963509ea914f3` |
| `rnor7` | AH109732 | `AHEnsDbs/v109/EnsDb.Rnorvegicus.v109.sqlite` | 179650560 | `2453c59e1054157798df5a4114650911d046b5f7795af91a6a2171614ba0289b` |
| `rnor7.113` | AH106786 | `ensembl/release-107/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit` | 661989722 | `73eacccd66a8843cced837d0cb0a272b42b7214b4a7c62a3c48963509ea914f3` |
| `rnor7.113` | AH119437 | `AHEnsDbs/v113/EnsDb.Rnorvegicus.v113.sqlite` | 182464512 | `ed63640b3d6d4cd7d26f1410a0746c341a17f7d8028df03db8dc772e77faff89` |

The other groups:

| Group | Key | URL | Size | SHA-256 |
|---|---|---|---|---|
| `example` | `sirna` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184929/suppl/GSE184929_Fig3-TTR-2_groupD1_30mkg.vs.Mock.DESeq2.txt.gz` | 1124255 | `c418c8aebd2138b892d91a0ca0a99952b083dd71f1043533d2c35045b2551142` |
| `example` | `sirna` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184929/suppl/GSE184929_Fig3-TTR-2_groupD4_30mkg.vs.Mock.DESeq2.txt.gz` | 1085045 | `a2e39fb9215aa647ffd8b6e7e4e491fffccf438e8f706e826e24c88313624696` |
| `example` | `sirna` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184929/suppl/GSE184929_Fig3-TTR-1_groupD1_10mkg.vs.Mock.DESeq2.txt.gz` | 1094976 | `97598daf5670d94875149cca0a3afebb7309758e41c99c35bb10274498186248` |
| `mirdb` | `6.0` | `https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz` | 59483982 | `446636d2182212f44ff91623f1429f7922f9575b82226962a6213d592a7e109d` |
| `trrust` | `human` | `https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv` | 297659 | `9b909319ccc8e36588b5a1bd3640e0dfba03936c6fc154cc4304d742bec11bed` |
| `trrust` | `mouse` | `https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv` | 220353 | `4847c960c7c2ce0173b056184137b1cc71e8d61d1b1edc825247504ad882dd55` |
| `string` | `9606` | `https://stringdb-downloads.org/download/protein.links.v12.5/9606.protein.links.v12.5.txt.gz` | 62813591 | `3d6894ef3a0a8578b70421b7d6bac159c41ed2f6c1b5da21ab7fabb4ac21cd25` |
| `string` | `10090` | `https://stringdb-downloads.org/download/protein.links.v12.5/10090.protein.links.v12.5.txt.gz` | 75044046 | `b3318ac1cbd0fe0e81a3d6b0b432352f72733b7cf36d2437423a5cbf2a98e9d0` |
| `string` | `10116` | `https://stringdb-downloads.org/download/protein.links.v12.5/10116.protein.links.v12.5.txt.gz` | 75736490 | `4d53bee2bc10daa046aae67532152c7389ba2374dacb18cf2b2cd3886b253387` |

The URLs of R: the GEO tables are at `R/get_data.R:23`, `:27`, and `:31`, miRDB at `R/get_data.R:39`, and TRRUST at `R/get_network_matches.R:21-22`. R reads STRING through `STRINGdb::STRINGdb$new` (`R/get_stringPPI_to_stringGGI.R:15`). Without a version, STRINGdb 2.22.0 uses the current version of the STRING API, 12.5 on 2026-11-20. The registry pins the version 12.5 (`D-041`).

## Review Focus

1. A connection that stops at 60% of a large file: no final file and no `.part` file, a `FetchError` that names the URL and the received byte count, and a second run that starts again. Owner: Task 2.
2. A `dest` that already holds the final file: no request when its digest agrees. When its digest differs, a `FetchError` that names the path and the two digests, and the file stays. Owner: Task 2.
3. A server that sends more bytes than the registry size, or no `Content-Length`: the read stops at the registry size plus 1 byte, and the chunk size limits the memory. Owner: Task 2.
4. A redirect from `https` to `http`: a `FetchError` before any request to the `http` URL. Owner: Task 2.
5. `seedmatch-fetch reference <build> DEST`, then `Reference.from_directory(DEST)`: the reference opens, because DEST has one `.sqlite` file and one `.2bit` file. Owner: Task 3.

---

### Task 1: The registry

**Files:**
- Make: `src/seedmatch/resources/__init__.py`, `src/seedmatch/resources/registry.py`
- Make: `tests/resources/test_registry.py`
- Modify: `tests/annotation/test_refs.py` (the registry gives the rnor7 file names and digests)

**Interfaces:**
- Produces, in `seedmatch.resources`:

```python
class Group(StrEnum):
    REFERENCE = "reference"
    EXAMPLE = "example"
    MIRDB = "mirdb"
    TRRUST = "trrust"
    STRING = "string"

@dataclass(frozen=True)
class Resource:
    url: str
    size: int
    sha256: str
    @property
    def name(self) -> str: ...  # the last part of the URL path

@dataclass(frozen=True)
class Entry:
    group: Group
    key: str
    description: str
    resources: tuple[Resource, ...]

REGISTRY: Final[tuple[Entry, ...]]
STRING_VERSION: Final = "12.5"

def lookup(group: Group, key: str, *, registry: Sequence[Entry] = REGISTRY) -> Entry: ...
```

- `lookup` raises `InputError` for an unknown key. The message names the group, the key, and the permitted keys of the group in registry order. A `group` that is not a `Group` member raises `InputError`, with `as_member` of `seedmatch._checks`.
- A reference entry has two resources: the 2bit file, then the EnsDb file. Its description names the AnnotationHub IDs, for example `AH106786 (2bit) and AH109732 (EnsDb), Ensembl 109`.

- [ ] **Step 1: Write the failing tests.**
  - `test_each_resource_has_an_https_url_a_size_and_a_digest`: the scheme is `https`, the size is above 0, and `sha256` is 64 lowercase hexadecimal digits.
  - `test_the_name_is_the_last_part_of_the_url_path`.
  - `test_each_group_and_key_occurs_one_time`, and `test_the_names_in_one_entry_are_different`.
  - `test_the_reference_keys_are_the_builds_of_r`: the keys of `Group.REFERENCE` are `hg38`, `mm39`, `mm10`, `rnor6`, `rnor7`, and `rnor7.113`, in this order. Each entry has a `.2bit` resource, then a `.sqlite` resource.
  - `test_lookup_names_the_permitted_keys`: `lookup(Group.TRRUST, "rat")` raises `InputError` with the exact message `trrust key 'rat' is not in the registry; expected one of 'human', 'mouse'`.
  - `test_lookup_rejects_text_for_the_group`.
  - `test_rnor7_and_rnor7_113_share_the_2bit_resource`.
  - In `tests/annotation/test_refs.py`, `ENSDB`, `TWOBIT`, and their rows of `RNOR7_FILES` come from `lookup(Group.REFERENCE, "rnor7")`. The GTF row stays, because the registry has no GTF. Thus `test_rnor7_files_are_the_pinned_files` makes sure that the registry values agree with the local files.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/resources -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.resources'`.
- [ ] **Step 3: Write `registry.py` and `__init__.py`.** Copy the values of "The registry" verbatim. `__init__.py` exports the names of the Interfaces block.
- [ ] **Step 4: Run the tests.** Run `uv run pytest tests/resources tests/test_package.py -q`, `uv run ruff check src tests`, `uv run ruff format --check src tests`, and `uv run mypy`. Expected: PASS. Then run `SEEDMATCH_TEST_REFS=~/Development/refs/seedmatch uv run pytest tests/annotation/test_refs.py -q -k pinned`. Expected: 3 passed, not skipped.
- [ ] **Step 5: Checkpoint.**

### Task 2: The download

**Files:**
- Make: `src/seedmatch/fetch/__init__.py`, `src/seedmatch/fetch/download.py`, `src/seedmatch/fetch/http.py`
- Modify: `src/seedmatch/errors.py` (`FetchError`), `pyproject.toml` (the per-file ignore `"tests/fetch/**" = ["TID251"]`, because these tests import `seedmatch.fetch` and `urllib.request`)
- Make: `tests/fetch/test_download.py`, `tests/fetch/test_http.py`, `tests/fetch/memory_client.py`

**Interfaces:**
- Consumes: `Resource`, `Entry` (Task 1).
- Produces, in `seedmatch.fetch`:

```python
class HttpResponse(Protocol):
    status: int
    def header(self, name: str) -> str | None: ...
    def read(self, size: int) -> bytes: ...
    def close(self) -> None: ...

@dataclass(frozen=True)
class FetchOpts:
    open_url: Callable[[str], HttpResponse]  # one GET, with no redirect
    chunk_size: int = 1 << 20
    max_redirects: int = 5

DEFAULT_FETCH_OPTS: Final[FetchOpts]  # open_url sends a GET with urllib.request

def fetch_resource(resource: Resource, dest: str | os.PathLike[str], *, opts: FetchOpts = DEFAULT_FETCH_OPTS) -> Path: ...
def fetch_entry(entry: Entry, dest: str | os.PathLike[str], *, opts: FetchOpts = DEFAULT_FETCH_OPTS) -> tuple[Path, ...]: ...
```

- `FetchError(SeedMatchError)` goes into `seedmatch/errors.py`.

Rules of `fetch_resource`:
- `dest` is a directory. If it is absent, make it with its parents. If it is a file, raise `InputError`.
- Let `final` be `dest / resource.name`. If `final` exists, read it in chunks. If its size and digest agree, return `final` with no request. Otherwise, raise `FetchError` that names `final`, the expected digest, and the actual digest, and tells the caller to remove the file. Do not change the file.
- A URL that is not `https`, as the start URL or as a redirect target, raises `FetchError` before the request to that URL.
- A status of 301, 302, 303, 307, or 308 goes to the `Location` header, resolved against the current URL. More than `max_redirects` redirects raise `FetchError`. A status other than 200 and the redirects raises `FetchError` with the status and the URL.
- Write the body to `dest / (resource.name + ".part")` in chunks of `chunk_size`, and update the SHA-256 with each chunk. Stop and raise `FetchError` when the received bytes pass `resource.size`. A `.part` file of an earlier run is replaced.
- At the end of the body, a size or a digest that differs raises `FetchError` with the URL, the expected value, and the actual value.
- On each error, and on an `OSError` or a `KeyboardInterrupt` during the transfer, delete the `.part` file, close the response, and raise again. An `OSError` becomes a `FetchError` that names the URL and the received byte count, with the `OSError` as its cause.
- When the size and the digest agree, `os.replace` the `.part` file to `final`, and return `final`.
- `fetch_entry` fetches each resource of the entry in order, and returns the paths.

Rules of the default client (`http.py`):
- `open_url` sends a GET with an opener of `urllib.request` whose redirect handler does not go to the target of a redirect. A 3xx `HTTPError` becomes a response with its status and its headers. Another `HTTPError` also becomes a response with its status, thus `fetch_resource` gives the message.
- The request has the header `User-Agent: seedmatch/<version>`, with the version of the package, and a timeout of 60 s for each socket operation.
- The opener of the default client is a parameter of a private function. Thus a test can give it an `https` handler with canned responses and no socket.

- [ ] **Step 1: Write `tests/fetch/memory_client.py`.** `MemoryClient(routes: Mapping[str, Route])` is an `open_url` function. A `Route` gives the status, the headers, the body, and an optional byte offset where `read` raises `ConnectionResetError`. The client records each requested URL in `requests`.
- [ ] **Step 2: Write the failing tests.**
  - `test_fetch_writes_the_file_and_leaves_no_part_file`.
  - `test_fetch_makes_the_dest_directory`, and `test_fetch_rejects_a_dest_that_is_a_file`.
  - `test_a_present_file_with_the_digest_causes_no_request`: `requests` is empty.
  - `test_a_present_file_with_another_digest_is_an_error_and_stays`: the message names the path and the two digests, and the bytes of the file do not change.
  - `test_a_digest_mismatch_leaves_no_file`, and `test_a_short_body_leaves_no_file`: the message names the URL, the expected value, and the actual value.
  - `test_a_long_body_stops_at_the_size_plus_one_byte`: a body of 10 times the size, with no `Content-Length`. The client records the bytes that `read` returned, and they are at most `size + chunk_size`.
  - `test_an_interrupted_transfer_leaves_no_file`: the route raises `ConnectionResetError` at 60% of the body. The error is `FetchError`, its cause is the `ConnectionResetError`, and the message gives the received byte count. A second call with a good route writes the file.
  - `test_a_stale_part_file_is_replaced`.
  - `test_an_http_url_is_an_error_before_a_request`: `requests` is empty.
  - `test_a_redirect_to_https_is_followed`, `test_a_relative_redirect_is_resolved`, `test_a_redirect_to_http_is_an_error_before_the_request`, and `test_too_many_redirects_is_an_error`.
  - `test_a_404_is_an_error_with_the_status`.
  - `test_fetch_entry_returns_the_paths_in_order`.
  - `tests/fetch/test_http.py` uses an opener whose `https` handler serves canned responses. For a redirect, `open_url` gives the status 302 and the `Location`, with no request to the target. For a 404, it gives the status 404. It sends the `User-Agent` header, and it reads the body.
- [ ] **Step 3: Run the tests.** Run `uv run pytest tests/fetch -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.fetch'`.
- [ ] **Step 4: Write `errors.py`, `http.py`, `download.py`, and `__init__.py`.**
- [ ] **Step 5: Run the tests.** Run `uv run pytest tests/fetch tests/test_package.py -q`, ruff, and mypy. Expected: PASS. `tests/test_package.py` still shows that `import seedmatch` loads no `seedmatch.fetch` module.
- [ ] **Step 6: Checkpoint.**

### Task 3: The command and the documents

**Files:**
- Make: `src/seedmatch/fetch/cli.py`, `src/seedmatch/fetch/__main__.py`, `tests/fetch/test_cli.py`
- Modify: `pyproject.toml` (`[project.scripts]`), `README.md`, `docs/r-map.md`, `docs/design.md` ("Fetch behavior")

**Interfaces:**
- Consumes: `REGISTRY`, `Group`, `lookup` (Task 1), `fetch_resource`, `FetchOpts`, `DEFAULT_FETCH_OPTS`, `FetchError` (Task 2).
- Produces: `main(argv: Sequence[str] | None = None, *, opts: FetchOpts = DEFAULT_FETCH_OPTS, registry: Sequence[Entry] = REGISTRY) -> int` in `seedmatch.fetch.cli`, and the script `seedmatch-fetch = "seedmatch.fetch.cli:main"`.

Rules of the command (`argparse`):
- The subcommands are `reference BUILD DEST`, `example NAME DEST`, `mirdb DEST`, `trrust ORGANISM DEST`, `string TAXON DEST`, and `list`. `mirdb` takes the only key of its group.
- A fetch calls `fetch_resource` for each resource of the entry. A file that is at `dest / resource.name` before the call is `present`. For another file, the command prints its name and size to stderr before the call. After each call, it prints the path and `present` or `downloaded` to stdout, separated by a tab.
- `list` prints one line for each resource: the group, the key, the file name, the size in bytes, and the description, separated by tabs.
- An `InputError` or a `FetchError` prints `seedmatch-fetch: <message>` to stderr and returns 1. A usage error of `argparse` exits with 2. Success returns 0.
- `python -m seedmatch.fetch` runs `main`.

- [ ] **Step 1: Write the failing tests.**
  - `test_fetch_reference_then_open_it`: a test registry has a `reference` entry `syn` whose resources serve the bytes of `tests/data/annotation/synthetic.2bit` and `synthetic.sqlite` through `MemoryClient`. `main(["reference", "syn", str(dest)], opts=..., registry=...)` returns 0, and `Reference.from_directory(dest)` opens. The test URL names end with `.2bit` and `.sqlite`.
  - `test_a_second_run_reports_present`.
  - `test_list_prints_each_resource`: one line for each resource of the default registry, with tab-separated fields.
  - `test_an_unknown_key_returns_1_and_names_the_keys`.
  - `test_a_fetch_error_returns_1_with_the_message`.
  - `test_a_usage_error_exits_with_2`.
  - `test_the_script_is_declared`: `pyproject.toml` has `seedmatch-fetch = "seedmatch.fetch.cli:main"`.
  - `test_python_m_runs_the_command`: `subprocess.run([sys.executable, "-m", "seedmatch.fetch", "list"])` returns 0. It makes no network connection.
- [ ] **Step 2: Run the tests.** Expected: FAIL with an `ImportError` for `seedmatch.fetch.cli`.
- [ ] **Step 3: Write `cli.py`, `__main__.py`, and the script line.**
- [ ] **Step 4: Run the tests.** Run `uv run pytest tests/fetch tests/resources tests/test_readme.py tests/test_package.py -q`, ruff, and mypy. Expected: PASS.
- [ ] **Step 5: Write the documents.**
  - `README.md`: a "Resources" part. The sandbox has no network access, thus a machine with network access runs `seedmatch-fetch`, and the sandbox reads the files from a path. Give the commands for a reference, the example tables, miRDB, TRRUST, and STRING, and then `Reference.from_directory(DEST)` and `read_deseq2_table`. Tell the reader to use one DEST for each build.
  - `docs/r-map.md`: the `load_species_anno_db` row names `seedmatch-fetch reference`. New rows for `get_example_data` (`seedmatch-fetch example sirna` and `seedmatch-fetch mirdb`, and `read_deseq2_table` for the tables) and for the TRRUST URLs of `get_network_matches` and the STRING download of `get_stringPPI_to_stringGGI`.
  - `docs/design.md`, "Fetch behavior": the rules of Task 2 for a present file with another digest, for the redirects, and for the size limit of the read.
  - `docs/design.md`, "Registry": the R line references of the URLs are stale. Use the lines of "The registry" in this plan.
- [ ] **Step 6: Checkpoint.**

### Task 4: Verification of the plan

- [ ] **Step 1: Run the checks one time.** Run `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, and `uv run pytest`. Expected: each passes.
- [ ] **Step 2: Fetch the small files for real (controller).** On the machine with network access, run `uv run seedmatch-fetch example sirna`, `trrust human`, and `string 10116`, each into a new directory under `tmp/claude/fetch-check/`. Expected: exit 0, and a second run prints `present` for each file. Then read each GEO table with `read_deseq2_table`.
- [ ] **Step 3: Fetch one reference file for real (controller).** Run `uv run seedmatch-fetch reference rnor7.113 tmp/claude/fetch-check/rnor7.113`, and compare the 2bit file with the local rnor7 2bit file (`cmp`). Then remove the directory. Checkpoint.
