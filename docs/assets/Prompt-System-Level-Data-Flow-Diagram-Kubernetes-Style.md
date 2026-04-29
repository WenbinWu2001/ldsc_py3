## Prompt: System-Level Data Flow Diagram (Kubernetes-Style) for a Software Package

### Objective

Create a **system-level architecture diagram** that visualizes:

1. **Data flow across the package**
2. **Layered architecture (clear boundaries)**
3. **Public vs private interfaces**
4. **Module-level responsibilities (coarse-grained)**

The diagram should follow **Kubernetes-style architecture conventions**, not DAG-style or low-level function graphs.

------

### 1. Core Structure

#### 1.1 Layout

- Use a **left → right flow** for data movement
- Use **vertical stacking for layers**
- Enforce a **strict layered hierarchy**

Example mental model:

```
[ External Input ]
        ↓
┌──────────────────────────────┐
│        CLI Layer (Public)    │
└──────────────────────────────┘
        ↓
┌──────────────────────────────┐
│ Workflow / API Layer (Public)│
└──────────────────────────────┘
        ↓
┌──────────────────────────────┐
│     Compute Kernel (Private) │
└──────────────────────────────┘
        ↓
[ Outputs / Files ]
```

------

### 2. Layer Representation (Critical)

#### 2.1 Each layer must be a **bounded container box**

- Use large rectangles to represent layers
- Label each layer clearly:

Format:

```
<Layer Name>
(public | private)
<module path>
```

Example:

```
CLI Dispatcher
(public)
ldsc.cli
```

------

#### 2.2 Enforce boundary semantics

- Arrows must **only flow downward across layers**
- No direct jumps skipping layers unless explicitly justified
- Private layers should **not expose outward arrows except through public layers**

------

### 3. Node Design (Inside Layers)

#### 3.1 Nodes represent **coarse modules or responsibilities**

NOT:

- individual functions
- implementation details

Use:

- "Parse CLI args"
- "Build config"
- "Run LD score computation"

------

#### 3.2 Node naming convention

Use:

```
<Verb> + <Object>
```

Examples:

- Load SNP data
- Intersect annotations
- Compute LD scores
- Write outputs

------

#### 3.3 Node styling

- Rectangles (uniform shape)
- Minimal text (≤ 3–5 words)
- No code inside nodes

------

### 4. Arrows (Data Flow)

#### 4.1 Arrow semantics

- Arrows represent **data movement**, not control logic
- Always directional (no bidirectional arrows unless necessary)

------

#### 4.2 Arrow labeling (optional but preferred)

Label important transitions with data types:

Examples:

- `VCF / SNP list`
- `BED annotations`
- `LD matrix`
- `.annot.gz`

------

### 5. Public vs Private Distinction

#### 5.1 Visual encoding

Use one of:

- Different border styles:
  - solid = public
  - dashed = private
- OR color coding:
  - darker = public
  - lighter = internal

------

#### 5.2 Explicit rule

- Only

   

  public layers

   

  may connect to:

  - external inputs
  - external outputs

------

### 6. External Systems / Inputs / Outputs

Represent external elements as **separate boxes outside layers**

Examples:

- Input data files
- User CLI invocation
- Output artifacts

Style:

- Simple rectangles
- Positioned at far left (inputs) and far right (outputs)

------

### 7. Grouping and Subsystems

If needed:

- Add **sub-boxes inside layers** for logical grouping
- Keep grouping shallow (max depth = 2)

Example:

```
Workflow Layer
 ├── Config handling
 ├── Job orchestration
 └── Output management
```

------

### 8. What to Avoid (Important Constraints)

Do NOT:

- Mix abstraction levels (functions + modules + systems)
- Draw full DAG graphs (too detailed)
- Include loops unless essential
- Overcrowd with text
- Show implementation details (algorithms, equations)

------

### 9. Level of Detail

Target:

- 8–20 nodes total
- 3–5 layers
- Readable in a single screen

This is a **system diagram**, not a code map.

------

### 10. Output Format (for rendering tools)

Specify that the diagram should be generated in one of:

- Mermaid (`flowchart LR` with subgraphs for layers)
- Graphviz (clustered subgraphs)
- Draw.io / Figma style layout

------

### 11. Example Input (to adapt)

Provide your layer table:

```
| Layer | Public? | Location |
| --- | --- | --- |
| CLI dispatcher | yes | ldsc.cli |
| workflow/config/output modules | yes | src/ldsc/*.py |
| compute kernel | no | src/ldsc/_kernel/*.py |
```

The model should:

- Convert this into **layer containers**
- Populate each layer with **representative nodes**
- Connect them via **data flow arrows**

------

### 12. User-Specific Additions For Reusable Workflow Documentation

Use this section when the diagram is meant to support package documentation rather than only produce a standalone picture.

#### 12.1 File-contract coverage

For each **major user-facing functionality**:

- Summarize the **required input files**
- Summarize the **main output files**
- Start from the **user-provided input files**
- End at the **program-emitted output files**

Where helpful, include a short example format for the files at the two ends of the flow.

Examples:

- input: ``CHR POS SNP CM annot1``
- input: ``SNP A1 A2 P N``
- output: ``.annot.gz``
- output: ``.l2.ldscore.gz``
- output: ``.h2.tsv``

#### 12.2 Streamlined per-feature flow

For each user-facing feature, describe a **streamlined flow**:

- keep only the main path
- omit edge cases unless they are architecturally essential
- show the skeleton of the workflow rather than every branch

The goal is to make it easy for a user to answer:

- what files do I need to provide?
- what layers/modules process them?
- what files do I get back?

#### 12.3 Layer and module traceability

Each flow should explicitly indicate:

- the **layer structure** used in that flow
- the main **modules** responsible for each stage

This should remain coarse-grained. Prefer module names such as:

- `ldsc.cli`
- `ldsc.ldscore_calculator`
- `ldsc.sumstats_munger`
- `ldsc.regression_runner`
- `ldsc._kernel.ldscore`

Avoid listing individual helper functions unless they are the public entry point.

#### 12.4 Preferred deliverable shape

When documenting a package with multiple workflows:

- Prefer **one overall package overview figure** covering all major functionalities
- Add **one or more supporting figures** only if the overview becomes too crowded
- Keep the figures **concise**
- Keep the figures **simple skeletons rather than dense diagrams**
- Do not overcrowd the figure with implementation detail or long labels

If a companion markdown document is produced alongside the figure, keep that document concise as well.

------

### 13. Expected Output Characteristics

The final diagram should:

- Look like a **Kubernetes architecture diagram adapted to software modules**
- Clearly show:
  - boundaries
  - interfaces
  - data movement
- Be understandable without reading the codebase
