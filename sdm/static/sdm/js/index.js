document.addEventListener('DOMContentLoaded', () => {
    $(window).scroll(function() {
        if ($(this).scrollTop() < 50) {
            $('.navbar').css('opacity', 0.0)
        } else {
            $('.navbar').css('opacity', 1.0)
        }
    });
});
