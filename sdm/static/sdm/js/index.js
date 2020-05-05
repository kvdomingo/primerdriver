var app_data = [{
    key: 0,
    name: 'Characterization',
    src: 'https://res.cloudinary.com/kdphotography-assets/image/upload/c_scale,w_500/v1/primerdriver/characterize.png',
    href: '#characterize',
    color: 'danger'
}, {
    key: 1,
    name: 'DNA-based',
    src: 'https://res.cloudinary.com/kdphotography-assets/image/upload/c_scale,w_500/v1/primerdriver/dna.png',
    href: '#dna-based',
    color: 'primary'
}, {
    key: 2,
    name: 'Protein-based',
    src: 'https://res.cloudinary.com/kdphotography-assets/image/upload/c_scale,w_500/v1/primerdriver/protein.png',
    href: '#protein-based',
    color: 'success'
}];

function parseEscapedJson(s) {
    return s.replace(/\\u[0-9a-fA-F]{4}/gi, function (match) {
        return String.fromCharCode(parseInt(match.replace(/\\u/g, ""), 16));
    });
}

document.addEventListener('DOMContentLoaded', function () {
    // $(window).scroll(function() {
    //     if ($(this).scrollTop() < 50) {
    //         $('.navbar').css('opacity', 0.0)
    //     } else {
    //         $('.navbar').css('opacity', 1.0)
    //     }
    // });

    expressionList = JSON.parse(parseEscapedJson($('#expression-systems').data('list')));
    parseEscapedJson($('#expression-systems').data('list'));

    ReactDOM.render(React.createElement(App, { stations: app_data, expressionList: expressionList.data }), document.getElementById('app'));
});