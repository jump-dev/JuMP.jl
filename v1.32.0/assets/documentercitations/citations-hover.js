/* Popups showing a bibliography entry when hovering over a citation link —
 * DocumenterCitations.
 *
 * The rendered entries are in `citations-data.js`, generated at build time by
 * the `WriteHoverData` pipeline step. It defines
 *
 *     window.DocumenterCitationsHoverData = {
 *       "<anchor name>": {"html": "<rendered entry>", "page": "<page path>"},
 *       ...
 *     }
 *
 * where `page` is the path of the page with the canonical bibliography block,
 * relative to the root of the site. A link is given a popup only if both its
 * fragment and its target page match an entry, so that links to something
 * other than a bibliography entry are left alone.
 *
 * The popup element is created here at runtime, and the entries come from the
 * data file rather than from the page. Hidden markup in the page (the usual
 * way to do this) would be invisible only to browsers that apply CSS: in a
 * text-mode browser, every page with citations would end up carrying a copy of
 * the bibliography.
 */
(function () {
  "use strict";

  // Must be read while the script is evaluated, not from a later callback
  var SCRIPT_SRC = document.currentScript && document.currentScript.src;

  var SHOW_DELAY = 150; // ms hovering the link before the popup appears
  var HIDE_DELAY = 300; // ms after leaving link and popup before it hides
  var POPUP_ID = "documentercitations-hover-popup";

  var popup = null; // singleton .citation-hover-popup element
  var showTimer = null;
  var hideTimer = null;
  var trigger = null; // the link the popup currently belongs to
  var pageURL = {}; // anchor name -> absolute URL of the bibliography page
  var pageKey = {}; // anchor name -> that URL in the form used for comparison
  var jumpedTo = null; // citation whose jump-focus was already skipped

  // A form of the URL that identifies the page, for comparing a link with the
  // recorded path of a bibliography page: without the fragment, and with a
  // trailing "/" resolved to "index.html", so that the comparison holds
  // irrespective of the `prettyurls` setting. This is not a valid URL (for a
  // `file://` page, `origin` is the string "null"), so it must not be used as
  // a base for resolving relative links.
  function pageKeyOf(url) {
    var pathname = url.pathname.endsWith("/") ? url.pathname + "index.html" : url.pathname;
    return url.origin + pathname + url.search;
  }

  function getPopup() {
    if (!popup) {
      popup = document.createElement("div");
      popup.className = "citation-hover-popup";
      popup.id = POPUP_ID;
      popup.setAttribute("role", "tooltip");
      // Keep the popup open while it is hovered or contains the focus, so that
      // it can be scrolled and links inside it can be clicked
      popup.addEventListener("mouseenter", cancelHide);
      popup.addEventListener("mouseleave", scheduleHide);
      popup.addEventListener("focusin", cancelHide);
      popup.addEventListener("focusout", scheduleHide);
      document.body.appendChild(popup);
    }
    return popup;
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

  // Best effort: the popup content is inserted after the page has been
  // typeset, so any math in the entry needs to be rendered again.
  function typesetMath(el, done) {
    if (window.MathJax && window.MathJax.typesetPromise) {
      window.MathJax.typesetPromise([el]).then(done);
    } else if (window.MathJax && window.MathJax.Hub) {
      window.MathJax.Hub.Queue(["Typeset", window.MathJax.Hub, el], done);
    } else if (
      typeof window.require === "function" &&
      window.require.specified &&
      window.require.specified("katex-auto-render")
    ) {
      window.require(["katex-auto-render"], function (renderMathInElement) {
        renderMathInElement(el, {
          delimiters: [
            { left: "$", right: "$", display: false },
            { left: "$$", right: "$$", display: true },
            { left: "\\[", right: "\\]", display: true },
          ],
        });
        done();
      });
    }
  }

  function showPopupFor(link) {
    var key = link.hash.slice(1);
    var entry = window.DocumenterCitationsHoverData[key];
    if (!entry) return;
    var p = getPopup();
    p.innerHTML = entry.html;
    // The entry was rendered relative to the bibliography page, so any link in
    // it must be resolved against that page, not the current one
    p.querySelectorAll("a[href]").forEach(function (a) {
      a.href = new URL(a.getAttribute("href"), pageURL[key]).href;
    });
    p.style.display = "block";
    trigger = link;
    link.setAttribute("aria-describedby", POPUP_ID);
    positionPopup(link);
    typesetMath(p, function () {
      if (trigger === link) positionPopup(link); // size may have changed
    });
  }

  function hidePopup() {
    if (trigger) {
      trigger.removeAttribute("aria-describedby");
      trigger = null;
    }
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

  // Keyboard focus shows the popup without delay, so that a screen reader
  // picks up the `aria-describedby` description when it announces the link.
  //
  // Following a link *to* a citation, e.g. one of the backlinks at the end of
  // a bibliography entry, also focuses it. That must not open a popup: the
  // reader asked to be taken to the citation, not to be shown the entry they
  // are coming from. Only the first focus is skipped, so that tabbing back to
  // the link later still describes it.
  function onLinkFocus(e) {
    var link = e.currentTarget;
    // The hash is up to date by the time the focus event fires (`hashchange`
    // only comes afterwards), so the link being the current fragment target
    // identifies the focus as the one that follows the jump
    if (link.id && "#" + link.id === window.location.hash && link.id !== jumpedTo) {
      jumpedTo = link.id;
      return;
    }
    cancelHide();
    if (showTimer) clearTimeout(showTimer);
    showPopupFor(link);
  }

  function init() {
    var data = window.DocumenterCitationsHoverData;
    if (!data || !SCRIPT_SRC) return;
    // The script is at <root>/assets/documentercitations/citations-hover.js
    var root = new URL("../..", SCRIPT_SRC);
    Object.keys(data).forEach(function (key) {
      var url = new URL(data[key].page, root);
      pageURL[key] = url.href;
      pageKey[key] = pageKeyOf(url);
    });
    var found = false;
    Array.prototype.forEach.call(document.links, function (link) {
      var key = link.hash.slice(1);
      if (!key || !data[key]) return;
      if (pageKeyOf(new URL(link.href)) !== pageKey[key]) return;
      link.addEventListener("mouseenter", onLinkEnter);
      link.addEventListener("mouseleave", onLinkLeave);
      link.addEventListener("focus", onLinkFocus);
      link.addEventListener("blur", onLinkLeave);
      found = true;
    });
    if (!found) return; // no citations on this page
    // The jump fires `hashchange` right after the focus it caused, so the
    // state that focus has just set must survive it: the reader is still at
    // that citation, and tabbing back to it should describe it. Moving on to
    // another fragment does clear the state, so that returning to the citation
    // later is again recognized as a jump.
    window.addEventListener("hashchange", function () {
      if ("#" + jumpedTo !== window.location.hash) jumpedTo = null;
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
