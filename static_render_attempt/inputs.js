// inputs.js
export const Inputs = {
    checkbox(options, config) {
        const container = document.createElement('div');
        container.className = 'input-container';
        
        if (config.label) {
            const label = document.createElement('h3');
            label.textContent = config.label;
            container.appendChild(label);
        }

        const form = document.createElement('form');
        
        options.forEach(option => {
            const label = document.createElement('label');
            const checkbox = document.createElement('input');
            checkbox.type = 'checkbox';
            checkbox.value = option;
            checkbox.checked = config.value.includes(option);
            
            label.appendChild(checkbox);
            label.appendChild(document.createTextNode(option));
            form.appendChild(label);
        });

        container.appendChild(form);
        
        // Return an object that mimics Observable's viewof structure
        return {
            element: container,
            value: config.value,
            addEventListener: (type, callback) => {
                form.addEventListener('change', (event) => {
                    const checkboxes = form.querySelectorAll('input[type="checkbox"]');
                    const selectedValues = Array.from(checkboxes)
                        .filter(cb => cb.checked)
                        .map(cb => cb.value);
                    this.value = selectedValues;
                    callback(selectedValues);
                });
            }
        };
    }
};
