function get_results(data_in) {
    const request = new XMLHttpRequest();
    request.open("POST", "/primer-design");
    request.onload = () => {
        const data = JSON.parse(request.responseText);
        data.forEach()
    }
}