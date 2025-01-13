export const initSvgExport = (containerId, filename = 'visualization') => {

    // Add export button to container
    const container = document.getElementById(containerId);
    const exportButton = document.createElement('button');
    exportButton.textContent = 'Export SVG';
    exportButton.className = 'export-button';
    container.appendChild(exportButton);
  
    const exportSvg = () => {
      // Get the SVG element
      const svgElement = container.querySelector('svg');
      if (!svgElement) {
        console.error('No SVG element found');
        return;
      }
  
      // Clone the SVG to avoid modifying the displayed version
      const clonedSvg = svgElement.cloneNode(true);
  
      // Ensure SVG has necessary attributes
      clonedSvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
      clonedSvg.setAttribute('version', '1.1');
  
      // Copy computed styles to inline styles
      const computedStyles = window.getComputedStyle(svgElement);
      clonedSvg.style.backgroundColor = computedStyles.backgroundColor;
  
      // Copy styles for all SVG elements
      const elements = clonedSvg.querySelectorAll('*');
      elements.forEach(el => {
        const computedStyle = window.getComputedStyle(el);
        const styleAttributes = [
          'fill',
          'stroke',
          'stroke-width',
          'opacity',
          'font-family',
          'font-size',
          'font-weight',
          'text-anchor',
          'dominant-baseline',
          'shape-rendering',
          'transform',
          'vector-effect'
        ];
        
        styleAttributes.forEach(attr => {
          if (computedStyle[attr]) {
            el.style[attr] = computedStyle[attr];
          }
        });
      });

      // Convert SVG to string
      const serializer = new XMLSerializer();
      const svgString = serializer.serializeToString(clonedSvg);
  
      // Create blob and download link
      const blob = new Blob([svgString], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = `${filename}.svg`;
      
      // Trigger download
      document.body.appendChild(link);
      link.click();
      
      // Cleanup
      document.body.removeChild(link);
      URL.revokeObjectURL(url);
    };
  
    // Add click handler
    exportButton.addEventListener('click', exportSvg);
  
    // Return cleanup function
    return () => {
      exportButton.removeEventListener('click', exportSvg);
      container.removeChild(exportButton);
    };
  };