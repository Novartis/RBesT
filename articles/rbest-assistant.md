# AI assistant skill

`RBesT` comes with an agent skill, `rbest-assistant`, that helps AI
assistants guide RBesT workflows: MAP priors from historical controls,
robustification, ESS, decision rules, operating characteristics and
probability of success. It asks before it codes and only uses verified,
runnable snippets.

**Please read before use.** The `RBesT` functions are tested and
validated; the AI agent’s output is not. The skill only guides the
agent: which code it writes, which inputs and assumptions it uses, and
how it interprets the results is decided by the agent alone. It is not
controlled by `RBesT`, it can differ between runs of the same prompt,
and it may be wrong. Treat it like a draft from a colleague: check the
code, the inputs and the conclusions, and rely only on results you have
reproduced by running the code yourself. Any analysis or design must be
reviewed by a qualified statistician. The assistant is experimental and
may change; watch for [updates](#updating), and send us your feedback
via [GitHub issues](https://github.com/Novartis/RBesT/issues).

## Getting started

Once the skill is installed (see below), try:

> Using the RBesT assistant skill: give me a short introduction to MAP
> priors. How can I use historical control data for a new trial? No
> analysis yet, just an overview. Then briefly tell me what else you can
> help me with in RBesT.

Naming the skill in plain words works with every assistant and install
type. To invoke it directly instead, type `/rbest:rbest-assistant` if
you installed the plugin (Claude Code, Claude Cowork, Posit Assistant,
VS Code), or `/rbest-assistant` if you installed the skill on its own
(in Codex: `$rbest-assistant`).

## Installation

Most routes install the skill straight from the GitHub repository; no
package release is needed:

- **Claude Cowork:** Customize → Plugins → Add marketplace →
  `Novartis/RBesT`, then install `rbest`.
- **Claude Code:** `/plugin marketplace add Novartis/RBesT`, then
  `/plugin install rbest@rbest`.
- **Claude Desktop and other apps that install skills by upload:**
  download
  [`rbest-assistant.zip`](https://github.com/Novartis/RBesT/releases/latest/download/rbest-assistant.zip)
  from the [latest GitHub
  release](https://github.com/Novartis/RBesT/releases/latest) and upload
  it as a skill (in Claude Desktop under Customize → Skills). Use this
  route if you cannot add a plugin marketplace. The zip holds the skill
  of that release, which can be older than the copy in the repository
  that the other routes install. Releases from `RBesT` 1.12-1 on have
  the zip attached.
- **ChatGPT Codex** (CLI, IDE extension or app): Codex reads skills from
  `~/.agents/skills` (for you) or `.agents/skills` (for one project).
  Run `npx skills add Novartis/RBesT`, use the R snippet under “Without
  plugins” below, or ask Codex’s built-in `$skill-installer` to install
  `inst/skills/rbest-assistant` from the `Novartis/RBesT` GitHub
  repository. Invoke it with `$rbest-assistant`.
- **Other agents:** `npx skills add Novartis/RBesT` (skill only).

Each installed version of `RBesT` also bundles a copy of the skill that
matches it, at
`system.file("skills/rbest-assistant", package = "RBesT")`.

**Any other tool: ask it to install the skill.** The skill uses the open
[Agent Skills](https://agentskills.io) format: a folder named
`rbest-assistant` that contains a `SKILL.md` file (YAML frontmatter with
`name` and `description`, then Markdown instructions) and a
`references/` folder of Markdown files. Most AI assistants that support
skills read this format. Paste this prompt into yours:

> Help me install the `rbest-assistant` agent skill. It is in the open
> Agent Skills format: a folder `rbest-assistant` containing `SKILL.md`
> (YAML frontmatter with `name` and `description`) and a `references/`
> folder. Its source is the GitHub repository `Novartis/RBesT`, folder
> `inst/skills/rbest-assistant`; it is also part of the installed RBesT
> R package at
> `system.file("skills/rbest-assistant", package = "RBesT")`. The
> repository root is also a Claude-format plugin marketplace
> (`.claude-plugin/marketplace.json`, plugin `rbest`). Tell me which
> install route you support (plugin, skill folder, or installer
> command), where you look for skills, walk me through installing it,
> and show me how to check that the skill is loaded. Keep the folder
> name `rbest-assistant`, since it must match the `name` in `SKILL.md`.

### Updating

- **Claude Code:** `/plugin marketplace update rbest`, or turn on
  auto-update for the marketplace under `/plugin` → Marketplaces.
- **Claude Cowork, Posit Assistant and VS Code:** use the update action
  of the app’s plugin settings (details below for Posit Assistant and VS
  Code).
- **Skill-only installs** (`npx skills`, Codex, a copied folder, an
  uploaded zip): repeat the install; for the zip, download it again from
  the [latest
  release](https://github.com/Novartis/RBesT/releases/latest) and
  replace the uploaded skill. The copy bundled with `RBesT` updates with
  the package.

### Positron, RStudio and VS Code

- **Positron / RStudio with Posit Assistant (recommended).** Posit
  Assistant reads the Claude marketplace format. Open ⋯ → Plugins (or
  `/plugin`) → Marketplaces, add `Novartis/RBesT`, then install `rbest`
  from Browse. You do not need `rbest-r`, because Posit Assistant runs R
  itself. To update, use “Check for updates” on the marketplace, then
  Update. The skill appears as `/rbest:rbest-assistant`. See the [Posit
  Assistant plugin
  docs](https://assistant.posit.co/docs/features/plugins/).

- **VS Code.** VS Code agent plugins read the Claude plugin and
  marketplace format. Add `Novartis/RBesT` as a marketplace (Chat: Open
  Customizations → Browse → Configure Marketplaces), or use the install
  links for
  [`rbest`](https://vscode.dev/redirect?url=vscode%3A%2F%2Fchat-plugin%2Finstall%3Fsource%3DNovartis%252FRBesT%26plugin%3Drbest)
  and
  [`rbest-r`](https://vscode.dev/redirect?url=vscode%3A%2F%2Fchat-plugin%2Finstall%3Fsource%3DNovartis%252FRBesT%26plugin%3Drbest-r)
  (deep link:
  `vscode://chat-plugin/install?source=Novartis%2FRBesT&plugin=rbest`,
  or `plugin=rbest-r`). This needs the `chat.plugins.enabled` setting.
  Updates come through “Check for and Apply Updates”. Enterprise Copilot
  policy may restrict marketplaces or MCP servers.

- **Without plugins**, in either editor, copy the skill bundled with
  your installed `RBesT` version:

  ``` r

  dir.create("~/.agents/skills", recursive = TRUE, showWarnings = FALSE)
  file.copy(system.file("skills/rbest-assistant", package = "RBesT"),
            "~/.agents/skills", recursive = TRUE)
  ```

  or run `npx skills add Novartis/RBesT`.

- **R connector without the plugin:** see btw’s [MCP
  documentation](https://posit-dev.github.io/btw/reference/mcp.html) for
  the `.vscode/mcp.json` or user `mcp.json` entry that runs
  `btw::btw_mcp_server()`.

## Running the code

**Results come from `RBesT` only.** The skill instructs the agent to
compute exclusively with the `RBesT` package in an R session, and never
to re-create its functionality in another language or with its own code.
What the agent delivers therefore depends on whether it can reach an R
with `RBesT`:

- **It can:** it runs the analysis there and reports the results
  together with the code.
- **It cannot:** it says so and offers either to guide you through the
  setup (R, `RBesT`, and a way for the agent to use them), or to write
  one self-contained R script that you download and run yourself. In
  that case it reports no results; you can bring the output back to have
  it interpreted.

If your organisation does not allow an agent to run R for you, or
requires analyses to run on a dedicated computing platform, tell the
agent. It then works in script-only mode from the start.

The skill names the few cases in which the agent may add numerical code
of its own in R on top of `RBesT`. At present there is one: assurance
for a design prior that links the two arms, which
[`pos2S()`](https://opensource.nibr.com/RBesT/reference/pos2S.md) cannot
express. It is obtained by quasi-Monte Carlo averaging of
[`oc2S()`](https://opensource.nibr.com/RBesT/reference/oc2S.md) results
(using a QMC package such as `spacefillr`) and is always reported with
its Monte Carlo standard error. Should an agent present results that it
computed in any other way, discard them and please [let us
know](https://github.com/Novartis/RBesT/issues).

**Recommended ways to have the code actually run:** Posit Assistant in
Positron or RStudio runs R in your live session, so it needs no
connector. Coding agents such as Claude Code or Copilot agent mode can
call `Rscript` themselves.

**Optional `rbest-r` connector.** This is mainly for Claude Cowork and
Claude Desktop, whose sandbox cannot reach your R. The `rbest-r` plugin
from the same marketplace starts btw’s MCP server so that the agent can
run R code in your local R installation. It needs a local R with
`RBesT`, `btw` and ideally `ragg`, and works on macOS and Linux. Install
`RBesT` as a ready-built binary rather than from source, preferably from
[Posit Public Package Manager](https://packagemanager.posit.co), which
also serves Linux binaries; the agent can guide you through it. It runs
R code without a sandbox; the app asks before each call. Organisations
may disable local MCP servers; the skill still works then, and the agent
hands over an R script instead. To get help with setup, ask your agent:

> Help me set up the optional `rbest-r` plugin from the `Novartis/RBesT`
> marketplace so my AI agent can run R code via the btw MCP server. I’m
> on \[macOS/Windows/Linux\]. Walk me through installing what’s needed,
> checking that the connector starts, and fixing it if Rscript isn’t
> found.

For Windows and manual setups, see btw’s [MCP
documentation](https://posit-dev.github.io/btw/reference/mcp.html).
