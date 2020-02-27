/* jshint esversion: 6 */

document.addEventListener('DOMContentLoaded', () => {

    $('#sequence').on('keyup keypress keydown mousedown mouseup change blur', function(e) {
        var cursorposition = $('#sequence').prop('selectionStart');
        document.getElementById('cursorposition').innerHTML = cursorposition + 1;
        document.getElementById('counter').innerHTML = `${this.value.length}/20000`;
    });

    $(document).on('keyup keypress keydown mousedown mouseup change blur', function(e) {
        var valid_sequence = ($('#sequence').val().length > 0);
        var valid_target = ($('#target').val().length > 0);
        var valid_position = !(isNaN(parseInt($('#position').val()))) && (parseInt($('#position').val()) > 0);
        var valid_replacement = ($('#target').val().length > 0);
        var valid_type = ($('#mut_type').val().length > 0);
        if (valid_sequence && valid_target && valid_position && valid_replacement && valid_type) {
            $('#submit').removeClass('disabled');
        }
        else {
            $('#submit').addClass('disabled');
        }
    });

    document.querySelector('#submit').onclick = () => {
        const xhttp = new XMLHttpRequest();
        const sequence = document.querySelector('#sequence').value;
        const target = document.querySelector('#target').value;
        const position = document.querySelector('#position').value;
        const replacement = document.querySelector('#replacement').value;
        const valued_settings = [
            '#Tm_range_min',
            '#Tm_range_max',
            '#gc_range_min',
            '#gc_range_max',
            '#length_min',
            '#length_max',
            '#flank5_range_min',
            '#flank5_range_max',
            '#flank3_range_min',
            '#flank3_range_max',
            '#primer_mode'
        ];
        const checked_settings = [
            '#terminate_gc',
            '#center_mutation'
        ];
        var settings = {};
        for (var set in valued_settings) {
            settings[set.slice(1)] = $(set).val();
        }
        for (set in checked_settings) {
            settings[set.slice(1)] = $(set)[0].checked;
        }
        settings = JSON.stringify(settings);

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
                <button class="btn btn-blue-grey btn-rounded" onclick="window.location.href='/dna-based'">Try again</button>
                `;
            document.querySelector('#form').innerHTML = contents;
        };
        const data = new FormData();
        data.append('mode', 'DNA');
        data.append('sequence', sequence);
        data.append('target', target);
        data.append('position', position);
        data.append('replacement', replacement);
        data.append('settings', settings);
        xhttp.send(data);
        return false;
    };
});
