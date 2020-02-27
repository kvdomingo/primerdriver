/* jshint esversion: 6 */

document.addEventListener('DOMContentLoaded', () => {
    document.getElementById('sequence').onkeyup = function () {
        document.getElementById('counter').innerHTML = `${this.value.length}/20000`;
    };

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
                <button class="btn btn-indigo btn-rounded" onclick="window.location.href='/characterize'">Back</button>
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
