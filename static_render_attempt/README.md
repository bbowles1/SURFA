# README

I made this page as a way to preview changes to the uORF rendering function in a way that's version-controlled.

To set up a local server via npm, I used:
1. `npm init` to initialize a project (which created package.json).
2. `npm install vite --save-dev` to install vite.
3. `npm run dev` to deploy a local webserver.


# Polish:

1. Add search bar to connect visual to python backend?
    - Need to display output in some sort of terminal if this is the route
2. Add filter for minimum uORF length.
3. Add checkbox for box labels (hard coded onto uORF boxes).
4. Park all scripts in a separate .js folder
5. Add xlabel "nucleotide position", add y-label for transcripts versus uORF names
8. Fix SVG render issues:
    - Purple CDS region is black in InkScape
9. Add export for uORF dataframe. 
10. Fix proper start codons for the python backend script.

# Stretch Goals:
1. Expand to other species / builds
2. Serve via API rather than Python backend
3. Convert to python package
4. Give all uORFs a unique gene-based identifier (ENSG.codon.num) to allow for comparison across transcripts
    - This would also allow annotation of other expressed transcripts, similar to a PEX measure
5. Export CSV of uORF annotations

# Demo Goals:
- Render an SVG and export it
- Develop a command line implementation that takes a transcript as the only arg
    - Meaning I must develop an importable database for all (human) transcript regions 
- CLI must have also have de novo implementation, which can generate a diagram based on a custom FASTA and GTF
    - The dataloader will ingest a FASTA and, if no matches are found, will raise a warning that the input FASTA seqids do not match
    - auto-format FASTA can attempt to format BED df IDs to match FASTA, but only if this is a simple change (ie appending chr prefix)
    - Warning that this may not handle MT sequences or edge cases like patches
- Must run on a Dockerfile
    - bedtools + python + npm

# Pipeline Steps
1. Build a Dockerfile with the following dependencies:
    - Bedtools
    - NPM
    - Pandas
    - UV (just for fun)
2. Use Dockerfile to serve .html to display a uorf.json file.
3. Use Python dataloader script to generate uORFs.json
4. Ensure Dockerfile can both render the uorfs.json for a HUMAN target, then call the html display module.

# Containerization
- Followed examples from https://github.com/astral-sh/uv-docker-example to set up the Dockerfile.
- Currently working with a uv Dockerfile in base directory.
- Can build image with `docker build . -t 'uorf-viewer'`.
- Can run with `docker run -it --rm --entrypoint /bin/bash -v /Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer:/mnt --name uorf_viewer a1b6bd367053`
- Copies all code to /app, python version is 3.13.
- I used uv init to begin a new project, `uv init uorf_viewer --bare`.
- Eventually, code should be placed in `src` dir.
- I added dependencies with `uv add pandas numpy`
- I added an install of npm, and ran `npm install` to download the dependencies from the package.json and package-lock.json files.
