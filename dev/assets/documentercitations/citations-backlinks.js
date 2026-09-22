/* Marking the backlink that leads back to the citation the reader came from —
 * DocumenterCitations.
 *
 * Every citation link carries the `id` that the backlinks at the end of the
 * corresponding bibliography entry point to. Following a citation therefore
 * identifies one particular backlink, but the browser cannot express that: the
 * link targets the entry as a whole, and nothing in the resulting page says
 * which of the entry's backlinks is the way back.
 *
 * So the id of the followed citation is remembered in `sessionStorage` (the
 * bibliography is usually on a different page, which rules out passing it in a
 * variable), and the backlink pointing to it is given the class
 * `is-current-citation`, styled in `citations.css`. Without JavaScript, or
 * when the reader arrives at the bibliography by other means, the backlinks
 * are simply all shown alike.
 */
(function () {
  "use strict";

  var CLASS = "is-current-citation";
  var STORAGE_KEY = "documentercitations-followed-citation";

  var marked = null; // the backlink currently carrying the class

  // `sessionStorage` throws in a browser that blocks storage, and is absent
  // for a page opened from the file system in some browsers
  function readFollowed() {
    try {
      return window.sessionStorage.getItem(STORAGE_KEY);
    } catch (e) {
      return null;
    }
  }

  function writeFollowed(id) {
    try {
      window.sessionStorage.setItem(STORAGE_KEY, id);
    } catch (e) {
      /* not being able to remember it only means no backlink is marked */
    }
  }

  function unmark() {
    if (marked) {
      marked.classList.remove(CLASS);
      marked = null;
    }
  }

  function mark(backlink) {
    unmark();
    if (!backlink) return;
    backlink.classList.add(CLASS);
    marked = backlink;
  }

  // The backlink pointing at the citation with the given id, if the page has
  // one. Any id that no backlink refers to gives `null`, which is what keeps
  // this from reacting to unrelated links.
  function findBacklink(id) {
    if (!id) return null;
    var backlinks = document.querySelectorAll(".citation-backlinks a");
    for (var i = 0; i < backlinks.length; i++) {
      if (backlinks[i].hash === "#" + id) return backlinks[i];
    }
    return null;
  }

  // The name of the anchor of the bibliography entry containing the given
  // backlink, i.e., the target of the citations of that reference
  function entryOf(backlink) {
    var element = backlink;
    while (element) {
      if (element.id) return element.id;
      element = element.parentNode;
    }
    return "";
  }

  // Whether following the given link means following a citation. Carrying an
  // `id` is not sufficient: `Documenter` gives one to some of its own
  // controls, e.g., the settings button, and clicking those must not be
  // remembered as a citation (nor unmark the backlink of an actual one). A
  // citation also links to the anchor of the bibliography entry.
  function isCitation(link) {
    return Boolean(link.id) && Boolean(link.hash);
  }

  // Mark the backlink for the citation remembered from an earlier click, but
  // only if the reader has in fact arrived at the entry that this citation
  // links to. Coming to the bibliography by other means — the sidebar, the
  // search, a link to a section — shows all backlinks alike, even though the
  // citation followed earlier in the session is still remembered.
  function markOnArrival() {
    var backlink = findBacklink(readFollowed());
    if (backlink && "#" + entryOf(backlink) !== window.location.hash) {
      backlink = null;
    }
    mark(backlink);
  }

  function init() {
    // A click on a citation records where to come back to
    document.addEventListener("click", function (e) {
      var element = e.target;
      while (element && element !== document) {
        if (element.tagName === "A") {
          if (isCitation(element)) {
            writeFollowed(element.id);
            // With the bibliography on the same page as the citation there is
            // no page load that would pick the id up again. The hash is only
            // updated after this handler has run, so `markOnArrival` cannot be
            // used here — but a click on a citation is arrival enough.
            mark(findBacklink(element.id));
          }
          return;
        }
        element = element.parentNode;
      }
    });
    markOnArrival();
    // Following a backlink, or any other jump within the page, leaves the mark
    // in place: it still shows where the reader came from. Coming back through
    // the history restores it, for a page served from the back/forward cache.
    window.addEventListener("pageshow", markOnArrival);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
