// main.js

import { render } from './render.js';
import { Inputs } from './inputs.js';
import UORFDetailsPanel from './display_details.js';
import { initSvgExport } from './exportSvg.js';

let currentVisualization = null;

// create the initial checkbox input
const selectCodons = Inputs.checkbox(
    ["ATG", "CTG", "TTG", "GTG", "TGG", "TCG", "UTU", "TTT", "TTC"],
    {
        label: "Start Codons",
        value: ["ATG", "CTG"]
    }
);

// Add the checkbox element to the controls
document.getElementById('controls').appendChild(selectCodons.element);

// func to load json from Flask API
function loadUorfData() {
    return fetch('/api/uorf').then(response => response.json());
  }


// func to calculate visualization dimensions
function getVisDimensions() {
    const visualizationDiv = document.getElementById('visualization');
    const width = Math.max(800, visualizationDiv.clientWidth);
    const height = Math.max(400, visualizationDiv.clientHeight);
    return { width, height };
}

// Function to update visualization based on selected codons
function updateVisualization(selectedCodons) {
    d3.select("#visualization").select("svg").remove();

    loadUorfData().then(uorfs => {
        const filteredRegions = Object.assign({}, uorfs, {
            regions: uorfs.regions.filter(region => 
                selectedCodons.includes(region.start_codon) || region.type === 'transcript'
            )
        });
        
        const { width, height } = getVisDimensions();
        
        const svg = d3.select("#visualization")
            .append("svg")
            .attr("width", width)
            .attr("height", height);
        
        currentVisualization = {
            data: filteredRegions,
            selectedCodons
        };
        
        render(svg.node(), width, height, filteredRegions, {
            onRegionHover: (region) => {
                ReactDOM.render(
                    React.createElement(UORFDetailsPanel, { region }),
                    document.getElementById('details-panel')
                );
            },
            onRegionLeave: () => {
                ReactDOM.render(
                    React.createElement(UORFDetailsPanel, { region: null }),
                    document.getElementById('details-panel')
                );
            }
        });
    }).catch(error => {
        console.error("Error loading the JSON file:", error);
    });
}

// Handle window resize
let resizeTimeout;
window.addEventListener('resize', () => {
    clearTimeout(resizeTimeout);
    resizeTimeout = setTimeout(() => {
        if (currentVisualization) {
            updateVisualization(currentVisualization.selectedCodons);
        }
    }, 250); // Debounce resize events
});

// Listen for changes to the checkboxes
selectCodons.addEventListener('change', updateVisualization);

// Initial render
updateVisualization(selectCodons.value);

// Initial render of the details panel
ReactDOM.render(
    React.createElement(UORFDetailsPanel, { region: null }),
    document.getElementById('details-panel')
);

// Initialize export functionality
const cleanup = initSvgExport('container', 'uorf-diagram');
