const CASE_LABELS = { pi: "PI Baseline", mco: "MCO 6-month", union: "PI + MCO" };
const STATUS_ORDER = ["done", "done-native-island", "processing", "partial", "none", "unknown"];
const STATUS_LABELS = {
  done: "Done",
  "done-native-island": "Native island",
  processing: "Processing",
  partial: "Partial",
  none: "None",
  unknown: "Unknown",
};

const elements = {
  caseButtons: [...document.querySelectorAll("[data-case]")],
  summaryCards: document.getElementById("summaryCards"),
  progressCaption: document.getElementById("progressCaption"),
  progressBar: document.getElementById("progressBar"),
  progressLegend: document.getElementById("progressLegend"),
  snapshotTime: document.getElementById("snapshotTime"),
  searchInput: document.getElementById("searchInput"),
  statusFilter: document.getElementById("statusFilter"),
  subtreeFilter: document.getElementById("subtreeFilter"),
  sortOrder: document.getElementById("sortOrder"),
  resultCount: document.getElementById("resultCount"),
  routineRows: document.getElementById("routineRows"),
  emptyState: document.getElementById("emptyState"),
};

let currentCase = "union";
let currentData = null;
const cache = new Map();

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>"']/g, character => ({
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#39;",
  })[character]);
}

function formatNumber(value) {
  return new Intl.NumberFormat("en-US").format(Number(value || 0));
}

function formatDate(value, includeTime = false) {
  if (!value) return "—";
  const options = includeTime
    ? { year: "numeric", month: "short", day: "numeric", hour: "2-digit", minute: "2-digit", timeZoneName: "short" }
    : { year: "numeric", month: "short", day: "numeric" };
  return new Intl.DateTimeFormat("en-US", options).format(new Date(value));
}

function percent(part, total) {
  return total ? ((part / total) * 100).toFixed(1) : "0.0";
}

async function loadCase(caseName) {
  setLoading(true);
  try {
    if (!cache.has(caseName)) {
      const response = await fetch(`data/${caseName}.json`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      cache.set(caseName, await response.json());
    }
    currentCase = caseName;
    currentData = cache.get(caseName);
    elements.caseButtons.forEach(button => button.classList.toggle("active", button.dataset.case === caseName));
    populateSubtrees();
    renderSummary();
    renderRows();
  } catch (error) {
    elements.summaryCards.innerHTML = `<div class="error-message">Unable to load the ${escapeHtml(CASE_LABELS[caseName])} snapshot: ${escapeHtml(error.message)}</div>`;
    elements.progressCaption.textContent = "Snapshot unavailable";
    elements.resultCount.textContent = "No data";
  } finally {
    setLoading(false);
  }
}

function setLoading(isLoading) {
  elements.caseButtons.forEach(button => { button.disabled = isLoading; });
  if (isLoading) elements.resultCount.textContent = "Loading…";
}

function renderSummary() {
  const { summary, generated_at: generatedAt } = currentData;
  const counts = summary.counts;
  const cards = [
    ["Total", formatNumber(summary.total), ""],
    ...STATUS_ORDER.map(status => [
      status,
      formatNumber(counts[status] || 0),
      `${formatNumber(summary.loc[status] || 0)} LOC`,
    ]),
  ];
  elements.summaryCards.innerHTML = cards.map(([label, value, detail]) => `
    <article class="summary-card">
      <span>${escapeHtml(label)}</span>
      <strong>${escapeHtml(value)}</strong>
      ${detail ? `<small>${escapeHtml(detail)}</small>` : ""}
    </article>
  `).join("");

  elements.progressCaption.textContent = `${CASE_LABELS[currentCase]} · ${formatNumber(summary.total)} routines`;
  elements.snapshotTime.textContent = `Exported ${formatDate(generatedAt, true)}`;
  elements.progressBar.innerHTML = STATUS_ORDER.map(status => {
    const count = counts[status] || 0;
    return `<span class="progress-segment status-${status}" style="width:${percent(count, summary.total)}%" title="${escapeHtml(STATUS_LABELS[status])}: ${formatNumber(count)}"></span>`;
  }).join("");
  elements.progressLegend.innerHTML = STATUS_ORDER.map(status => `
    <span class="legend-item">
      <span class="legend-dot status-${status}"></span>
      ${escapeHtml(STATUS_LABELS[status])} ${formatNumber(counts[status] || 0)}
    </span>
  `).join("");
}

function populateSubtrees() {
  const previous = elements.subtreeFilter.value;
  const subtrees = [...new Set(currentData.routines.map(row => row.subtree).filter(Boolean))].sort();
  elements.subtreeFilter.innerHTML = '<option value="">All subtrees</option>' + subtrees
    .map(subtree => `<option value="${escapeHtml(subtree)}">${escapeHtml(subtree)}</option>`)
    .join("");
  if (subtrees.includes(previous)) elements.subtreeFilter.value = previous;
}

function filteredRows() {
  const query = elements.searchInput.value.trim().toLowerCase();
  const status = elements.statusFilter.value;
  const subtree = elements.subtreeFilter.value;
  const rows = currentData.routines.filter(row => {
    const matchesQuery = !query || `${row.routine} ${row.relpath}`.toLowerCase().includes(query);
    return matchesQuery && (!status || row.display_status === status) && (!subtree || row.subtree === subtree);
  });

  const nullLast = value => value == null ? -1 : Number(value);
  const sorters = {
    path: (a, b) => `${a.relpath}:${a.routine}`.localeCompare(`${b.relpath}:${b.routine}`),
    status: (a, b) => STATUS_ORDER.indexOf(a.display_status) - STATUS_ORDER.indexOf(b.display_status) || a.routine.localeCompare(b.routine),
    "coverage-desc": (a, b) => nullLast(b.coverage_percent) - nullLast(a.coverage_percent),
    "codon-desc": (a, b) => nullLast(b.codon_coverage_percent) - nullLast(a.codon_coverage_percent),
    "updated-desc": (a, b) => String(b.updated_at || "").localeCompare(String(a.updated_at || "")),
  };
  return rows.sort(sorters[elements.sortOrder.value] || sorters.path);
}

function coverageCell(label, percentValue) {
  if (percentValue == null || label === "unknown") return "—";
  return escapeHtml(label || `${Number(percentValue).toFixed(1)}%`);
}

function renderRows() {
  if (!currentData) return;
  const rows = filteredRows();
  elements.resultCount.textContent = `${formatNumber(rows.length)} of ${formatNumber(currentData.routines.length)} routines`;
  elements.emptyState.hidden = rows.length !== 0;
  elements.routineRows.innerHTML = rows.map(row => {
    const status = STATUS_ORDER.includes(row.display_status) ? row.display_status : "unknown";
    return `
      <tr>
        <td class="routine-name">${escapeHtml(row.routine)}</td>
        <td class="source-path">${escapeHtml(row.relpath)}</td>
        <td><span class="status-pill status-${status}">${escapeHtml(STATUS_LABELS[status])}</span></td>
        <td class="numeric">${coverageCell(row.coverage_label, row.coverage_percent)}</td>
        <td class="numeric">${coverageCell(row.codon_coverage_label, row.codon_coverage_percent)}</td>
        <td class="date-cell">${escapeHtml(formatDate(row.updated_at))}</td>
      </tr>
    `;
  }).join("");
}

elements.caseButtons.forEach(button => button.addEventListener("click", () => loadCase(button.dataset.case)));
[elements.searchInput, elements.statusFilter, elements.subtreeFilter, elements.sortOrder]
  .forEach(control => control.addEventListener(control.tagName === "INPUT" ? "input" : "change", renderRows));

loadCase("union");
