
// Function to toggle source code visibility
document.addEventListener('DOMContentLoaded', function() {
    // Add toggle buttons for implementation sections
    const sections = document.querySelectorAll('.breathe-sectiondef');
    sections.forEach(function(section) {
        if (section.querySelector('.cpp-source')) {
            const btn = document.createElement('button');
            btn.textContent = 'Toggle Implementation';
            btn.className = 'toggle-impl-btn';
            btn.style.cssText = 'background: #2980b9; color: white; border: none; padding: 5px 10px; margin: 5px 0; cursor: pointer; border-radius: 3px;';
            btn.onclick = function() {
                const sources = section.querySelectorAll('.cpp-source');
                sources.forEach(function(src) {
                    src.style.display = src.style.display === 'none' ? 'block' : 'none';
                });
            };
            section.insertBefore(btn, section.firstChild);
        }
    });
    
    // Add a link to source browser in the navigation
    const nav = document.querySelector('.wy-nav-side .wy-menu-vertical');
    if (nav) {
        const sourceLi = document.createElement('li');
        sourceLi.className = 'toctree-l1';
        const sourceLink = document.createElement('a');
        sourceLink.href = '/source_browser.html';
        sourceLink.textContent = 'Source Code Browser';
        sourceLi.appendChild(sourceLink);
        nav.appendChild(sourceLi);
    }
});
