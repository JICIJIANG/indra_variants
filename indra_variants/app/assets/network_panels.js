/* Collapse/resize behavior for the network detail page's title header and
 * variant-map panel. Pure DOM manipulation via event delegation on
 * `document` so it keeps working across Dash's client-side page swaps and
 * needs no Dash callbacks/stores for what is purely visual state. */
(function () {
    function toggleHeader(btn) {
        var header = btn.closest('.net-header');
        if (!header) return;
        var body = header.querySelector('.net-header-body');
        if (!body) return;
        var handle = header.nextElementSibling;
        var collapsed = body.classList.toggle('is-collapsed');
        if (collapsed) {
            header.dataset.expandedFlex = header.style.flex || '';
            header.style.flex = '0 0 auto';
        } else {
            header.style.flex = header.dataset.expandedFlex || '';
        }
        if (handle && handle.classList.contains('header-resize-handle')) {
            handle.style.display = collapsed ? 'none' : '';
        }
        btn.textContent = collapsed ? '▸' : '▾';
    }

    function toggleMap(btn) {
        var panel = btn.closest('.net-map-panel');
        if (!panel) return;
        var body = panel.querySelector('.net-map-body');
        if (!body) return;
        var handle = panel.previousElementSibling;
        var collapsed = body.classList.toggle('is-collapsed');
        if (collapsed) {
            panel.dataset.expandedFlex = panel.style.flex || '';
            panel.classList.add('is-collapsed');
        } else {
            panel.classList.remove('is-collapsed');
            panel.style.flex = panel.dataset.expandedFlex || '';
        }
        if (handle && handle.classList.contains('network-resize-handle')) {
            handle.style.display = collapsed ? 'none' : '';
        }
        btn.textContent = collapsed ? '▸' : '▾';
    }

    function toggleSidebar(btn) {
        var sidebar = document.querySelector('.net-sidebar');
        var content = document.querySelector('.net-main-content');
        if (!sidebar) return;
        var collapsed = sidebar.classList.toggle('is-collapsed');
        if (content) content.classList.toggle('is-sidebar-collapsed', collapsed);
        btn.classList.toggle('is-collapsed', collapsed);
        btn.textContent = collapsed ? '▸' : '◂';
    }

    document.addEventListener('click', function (e) {
        var headerToggle = e.target.closest('.net-header-toggle');
        if (headerToggle) {
            toggleHeader(headerToggle);
            return;
        }
        var mapToggle = e.target.closest('.net-map-toggle');
        if (mapToggle) {
            toggleMap(mapToggle);
            return;
        }
        var sidebarToggle = e.target.closest('.sidebar-toggle');
        if (sidebarToggle) {
            toggleSidebar(sidebarToggle);
        }
    });

    document.addEventListener('mousedown', function (e) {
        var handle = e.target.closest('.network-resize-handle');
        if (!handle) return;
        var panel = handle.nextElementSibling;
        var container = handle.parentElement;
        if (!panel || !container) return;

        e.preventDefault();
        var prevCursor = document.body.style.cursor;
        document.body.style.cursor = 'row-resize';

        function onMove(ev) {
            var rect = container.getBoundingClientRect();
            var height = rect.bottom - ev.clientY;
            var min = 80;
            var max = Math.max(min, rect.height - 160);
            height = Math.max(min, Math.min(max, height));
            panel.style.flex = '0 0 ' + height + 'px';
        }

        function onUp() {
            document.removeEventListener('mousemove', onMove);
            document.removeEventListener('mouseup', onUp);
            document.body.style.cursor = prevCursor;
        }

        document.addEventListener('mousemove', onMove);
        document.addEventListener('mouseup', onUp);
    });

    document.addEventListener('mousedown', function (e) {
        var handle = e.target.closest('.header-resize-handle');
        if (!handle) return;
        var header = handle.previousElementSibling;
        if (!header || !header.classList.contains('net-header')) return;
        var titleEl = header.querySelector('.net-header-title');
        var subtitleEl = header.querySelector('.net-header-subtitle');
        if (!titleEl) return;

        e.preventDefault();
        var prevCursor = document.body.style.cursor;
        document.body.style.cursor = 'row-resize';

        var baseHeight = header.getBoundingClientRect().height;
        var headerStyle = getComputedStyle(header);
        var basePaddingTop = parseFloat(headerStyle.paddingTop);
        var basePaddingBottom = parseFloat(headerStyle.paddingBottom);
        var titleStyle = getComputedStyle(titleEl);
        var baseTitleSize = parseFloat(titleStyle.fontSize);
        var baseTitleMarginTop = parseFloat(titleStyle.marginTop);
        var baseTitleMarginBottom = parseFloat(titleStyle.marginBottom);
        var subtitleStyle = subtitleEl ? getComputedStyle(subtitleEl) : null;
        var baseSubtitleSize = subtitleStyle ? parseFloat(subtitleStyle.fontSize) : null;
        var baseSubtitleMarginBottom = subtitleStyle ? parseFloat(subtitleStyle.marginBottom) : null;
        var startY = e.clientY;

        function onMove(ev) {
            var newHeight = baseHeight + (ev.clientY - startY);
            var min = 32;
            var max = 220;
            newHeight = Math.max(min, Math.min(max, newHeight));
            var scale = Math.max(0.4, Math.min(1.8, newHeight / baseHeight));
            header.style.flex = '0 0 ' + newHeight + 'px';
            header.style.paddingTop = (basePaddingTop * scale) + 'px';
            header.style.paddingBottom = (basePaddingBottom * scale) + 'px';
            titleEl.style.fontSize = (baseTitleSize * scale) + 'px';
            titleEl.style.marginTop = (baseTitleMarginTop * scale) + 'px';
            titleEl.style.marginBottom = (baseTitleMarginBottom * scale) + 'px';
            if (subtitleEl) {
                subtitleEl.style.fontSize = (baseSubtitleSize * scale) + 'px';
                subtitleEl.style.marginBottom = (baseSubtitleMarginBottom * scale) + 'px';
            }
        }

        function onUp() {
            document.removeEventListener('mousemove', onMove);
            document.removeEventListener('mouseup', onUp);
            document.body.style.cursor = prevCursor;
        }

        document.addEventListener('mousemove', onMove);
        document.addEventListener('mouseup', onUp);
    });
})();
