# README

I made this page as a way to preview changes to the uORF rendering function in a way that's version-controlled.

To set up a local server via npm, I used:
1. `npm init` to initialize a project (which created package.json).
2. `npm install vite --save-dev` to install vite.
3. `npm run dev` to deploy a local webserver.


# Demo Goals:

1. Add search bar to connect visual to python backend.
2. Add filter for minimum uORF length.
3. Add checkbox for box labels (hard coded onto uORF boxes).
4. Add color-coding for start codon identity.
5. Add xlabel "nucleotide position", add y-label for transcripts versus uORF names
6. Adjust legend
    Opacity for CDS must match.
    Color coding must match transcript regions.
7. Add SVG export support.
8. Fix padding - it's cutting off top of graph
9. Add export for uORF dataframe. 
10. Fix proper start codons for the python backend script.

# Stretch Goals:
1. Expand to other species / builds
2. Serve via API rather than Python backend