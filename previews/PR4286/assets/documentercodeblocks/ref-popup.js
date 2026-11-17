/* Doxygen-style hover tooltips for reference links — DocumenterCodeBlocks.
 *
 * The build step embeds, per page, one hidden `.ref-tip` div per unique
 * reference target (signature + first-sentence brief) inside a
 * `.ref-tips` container — the same trick as doxygen's `.ttc` divs. Hovering
 * any `a.julia-ref` shows the tip for its href.
 *
 * Ambiguous references (several documented methods match; `data-ref-targets`
 * holds the candidate [label, href] pairs) show a clickable signature list
 * instead, with the hovered candidate's tip below it.
 *
 * Everything is static markup + this script; nothing is fetched.
 */
(function () {
  "use strict";

  var SHOW_DELAY = 150; // ms hovering the link before the tooltip appears
  var HIDE_DELAY = 300; // ms after leaving link+tooltip before it hides

  var popup = null; // singleton .ref-popup element
  var showTimer = null;
  var hideTimer = null;
  var tipMap = null; // href -> .ref-tip element

  function buildTipMap() {
    tipMap = new Map();
    document.querySelectorAll(".ref-tips .ref-tip").forEach(function (tip) {
      tipMap.set(tip.getAttribute("data-for"), tip);
    });
  }

  function getPopup() {
    if (!popup) {
      popup = document.createElement("div");
      popup.className = "ref-popup";
      popup.addEventListener("mouseenter", cancelHide);
      popup.addEventListener("mouseleave", scheduleHide);
      document.body.appendChild(popup);
    }
    return popup;
  }

  function parseTargets(link) {
    var raw = link.getAttribute("data-ref-targets");
    if (!raw) return null;
    try {
      var t = JSON.parse(raw);
      return Array.isArray(t) && t.length > 1 ? t : null;
    } catch (e) {
      return null;
    }
  }

  function tipClone(href) {
    var tip = tipMap.get(href);
    if (!tip) return null;
    var clone = tip.cloneNode(true);
    clone.removeAttribute("data-for");
    return clone;
  }

  function positionPopup(link) {
    var p = getPopup();
    var r = link.getBoundingClientRect();
    var margin = 8;
    var pw = p.offsetWidth;
    var ph = p.offsetHeight;
    var left = Math.max(margin, Math.min(r.left, window.innerWidth - pw - margin));
    var top = r.bottom + 6;
    if (top + ph > window.innerHeight - margin && r.top - 6 - ph > margin) {
      top = r.top - 6 - ph; // flip above when out of room below
    }
    p.style.left = left + window.scrollX + "px";
    p.style.top = top + window.scrollY + "px";
  }

  function showPopupFor(link) {
    var targets = parseTargets(link);
    var p = getPopup();
    p.innerHTML = "";

    if (targets) {
      // Disambiguation: candidate signature list + hovered candidate's tip.
      var list = document.createElement("ul");
      list.className = "ref-popup-targets";
      var detail = document.createElement("div");
      detail.className = "ref-popup-detail";
      targets.forEach(function (t, i) {
        var li = document.createElement("li");
        var a = document.createElement("a");
        a.textContent = t[0];
        a.href = new URL(t[1], location.href).href;
        li.appendChild(a);
        li.addEventListener("mouseenter", function () {
          selectTarget(list, li, detail, t[1], link);
        });
        list.appendChild(li);
        if (i === 0) selectTarget(list, li, detail, t[1], link);
      });
      p.appendChild(list);
      p.appendChild(detail);
    } else {
      // An arity-narrowed aggregate points at its variant tip via data-ref-tip
      // (several call sites can share the href but need different tips).
      var clone = tipClone(link.getAttribute("data-ref-tip") || link.getAttribute("href"));
      if (!clone) return; // no tip for this target — no tooltip
      p.appendChild(clone);
    }
    p.style.display = "block";
    positionPopup(link);
  }

  function selectTarget(list, li, detail, href, link) {
    list.querySelectorAll("li").forEach(function (x) {
      x.classList.remove("selected");
    });
    li.classList.add("selected");
    detail.innerHTML = "";
    var clone = tipClone(href);
    if (clone) {
      detail.appendChild(clone);
      positionPopup(link);
    }
  }

  function hidePopup() {
    if (popup) {
      popup.style.display = "none";
      popup.innerHTML = "";
    }
  }

  function cancelHide() {
    if (hideTimer) {
      clearTimeout(hideTimer);
      hideTimer = null;
    }
  }

  function scheduleHide() {
    cancelHide();
    hideTimer = setTimeout(hidePopup, HIDE_DELAY);
  }

  function onLinkEnter(e) {
    var link = e.currentTarget;
    cancelHide();
    if (showTimer) clearTimeout(showTimer);
    showTimer = setTimeout(function () {
      showPopupFor(link);
    }, SHOW_DELAY);
  }

  function onLinkLeave() {
    if (showTimer) {
      clearTimeout(showTimer);
      showTimer = null;
    }
    scheduleHide();
  }

  function init() {
    buildTipMap();
    if (tipMap.size === 0) return; // page has no reference targets
    document.querySelectorAll("a.julia-ref").forEach(function (link) {
      link.addEventListener("mouseenter", onLinkEnter);
      link.addEventListener("mouseleave", onLinkLeave);
      link.addEventListener("focus", onLinkEnter);
      link.addEventListener("blur", onLinkLeave);
    });
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape") hidePopup();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
