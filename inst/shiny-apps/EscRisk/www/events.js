$(function() {

    var editing_skeleton = 0;

    const skel_probs = document.getElementById('crm-skeleton');
    skel_probs.addEventListener('focusout', e => {
	editing_skeleton = editing_skeleton - 1;
	setTimeout(function() {
	    if (!editing_skeleton) {
		Shiny.setInputValue("editing_skeleton", false);
	    }
	}); // TODO: Any need to set non-zero delay?
    });
    skel_probs.addEventListener('focusin', e => {
	Shiny.setInputValue("editing_skeleton", true);
	editing_skeleton = editing_skeleton + 1;
    });

    const outputJ = document.getElementById('J');
    Shiny.addCustomMessageHandler('watch-J', function(J) {
        outputJ.innerHTML = J;
    });
});
