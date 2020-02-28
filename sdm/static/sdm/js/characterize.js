/* jshint esversion: 6 */

document.addEventListener('DOMContentLoaded', () => {

    $('#sequence').on('keyup keypress keydown mousedown mouseup change blur', function(e) {
        var cursorposition = $('#sequence').prop('selectionStart');
        document.getElementById('cursorposition').innerHTML = cursorposition + 1;
        document.getElementById('counter').innerHTML = `${this.value.length}/20000`;
    });

    $(document).on('keyup keypress keydown mousedown mouseup change blur', function(e) {
        var valid_sequence = ($('#sequence').val().length > 0);
        var valid_mismatch = ($('#mismatch').val().length > 0);
        var valid_type = ($('#mut_type').val().length > 0);
        if (valid_sequence && valid_mismatch && valid_type) {
            $('#submit').removeClass('disabled');
        }
        else {
            $('#submit').addClass('disabled');
        }
    });

    document.querySelector('#submit').onclick = () => {
        const xhttp = new XMLHttpRequest();
        const sequence = document.querySelector('#sequence').value;
        const mismatch = document.querySelector('#mismatch').value;
        const mut_type = document.querySelector('#mut_type').value;
        
        xhttp.open('POST', '/result');
        xhttp.onload = () => {
            const raw = xhttp.responseText;
            const data = JSON.parse(raw);
            var contents = `<div class="table-responsive text-nowrap"><table class="table table-sm table-bordered table-hover"><tbody>`;
            for (var key in data) {
                contents += `
                    <tr>
                        <th scope="row">${key}</th>
                        <td>${data[key]["1"]}</td>
                    </tr>`;
            }
            contents += `
                </tbody></table></div>
                <button class="btn btn-blue-grey btn-rounded" onclick="window.location.href='/characterize'">Try again</button>
                `;
            document.querySelector('#form').innerHTML = contents;
        };
        const data = new FormData();
        data.append('mode', 'CHAR');
        data.append('sequence', sequence);
        data.append('mismatched_bases', mismatch);
        data.append('mutation_type', mut_type);
        xhttp.send(data);
        return false;
    };
});
