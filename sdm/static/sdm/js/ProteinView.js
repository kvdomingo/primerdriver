var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var ProteinView = function (_React$Component) {
    _inherits(ProteinView, _React$Component);

    function ProteinView(props) {
        _classCallCheck(this, ProteinView);

        var _this = _possibleConstructorReturn(this, (ProteinView.__proto__ || Object.getPrototypeOf(ProteinView)).call(this, props));

        _this.genericNucleobaseHandler = function (e) {
            name = e.target.name;
            old_value = e.target.value.toUpperCase().split("");
            amino_acids = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
            value = [];
            old_value.forEach(function (char, i) {
                if (amino_acids.includes(char)) value.push(char);
            });
            value = value.join('');

            if (name === 'sequence') {
                var _this$setState;

                _this.setState((_this$setState = {}, _defineProperty(_this$setState, name, value), _defineProperty(_this$setState, 'cursorPosition', e.target.selectionStart + 1), _defineProperty(_this$setState, 'sequenceLength', e.target.value.length), _this$setState));
            } else {
                _this.setState(_defineProperty({}, name, value));
            }
        };

        _this.mutationTypeHandler = function (e) {
            if (e.target.value !== _this.state.mutation_type) {
                _this.setState({
                    target: '',
                    replacement: ''
                });
            }
            _this.setState({ mutation_type: e.target.value });
        };

        _this.genericSelectHandler = function (e) {
            name = e.target.name;
            value = e.target.value;
            _this.setState(_defineProperty({}, name, value));
        };

        _this.genericCheckedHandler = function (e) {
            name = e.target.name;
            value = e.target.checked;
            _this.setState(_defineProperty({}, name, value));
        };

        _this.genericTextHandler = function (e) {
            name = e.target.name;
            value = e.target.value;
            _this.setState(_defineProperty({}, name, value));
        };

        _this.genericFloatHandler = function (e) {
            name = e.target.name;
            value = parseFloat(e.target.value);

            if (name === 'Tm_range_min' && value >= _this.state.Tm_range_max) return;
            if (name === 'Tm_range_max' && value <= _this.state.Tm_range_min) return;
            if (name === 'gc_range_min' && value >= _this.state.gc_range_max) return;
            if (name === 'gc_range_max' && value <= _this.state.gc_range_min) return;

            return _this.setState(_defineProperty({}, name, value));
        };

        _this.genericIntHandler = function (e) {
            name = e.target.name;
            value = parseInt(e.target.value);
            _this.setState(_defineProperty({}, name, value));
        };

        _this.resetForm = function () {
            _this.setState(Object.assign({}, _this.formResetDefaults));
        };

        _this.formValidator = function () {
            validSequence = _this.state.sequenceLength > 0;
            validMutation = _this.state.mutation_type !== '';
            validMutationCode = _this.state.target.length > 0 || _this.state.replacement.length > 0;
            validSequence && validMutation && validMutationCode ? _this.setState({ submitValid: true }) : _this.setState({ submitValid: false });
            _this.formSettings = {
                Tm_range_min: _this.state.Tm_range_min,
                Tm_range_max: _this.state.Tm_range_max,
                gc_range_min: _this.state.gc_range_min,
                gc_range_max: _this.state.gc_range_max,
                length_min: _this.state.length_min,
                length_max: _this.state.length_max,
                flank5_range_min: _this.state.flank5_range_min,
                flank5_range_max: _this.state.flank5_range_max,
                flank3_range_min: _this.state.flank3_range_min,
                flank3_range_max: _this.state.flank3_range_max,
                forward_overlap5: _this.state.forward_overlap5,
                forward_overlap3: _this.state.forward_overlap3,
                terminate_gc: _this.state.terminate_gc,
                center_mutation: _this.state.center_mutation,
                primer_mode: _this.state.primer_mode,
                expression_system: _this.state.expression_system
            };
        };

        _this.submitHandler = function (e) {
            _this.setState({ loading: true });
            xhttp = new XMLHttpRequest();
            xhttp.open('POST', '/api');
            xhttp.onload = function () {
                if (xhttp.status === 200) {
                    res = JSON.parse(xhttp.responseText);
                } else {
                    res = 'Request failed. Please try again later.';
                    console.log(xhttp.statusText);
                    console.log(xhttp.response);
                    console.log(xhttp.responseText);
                }
                _this.props.responseCatcher(e, res, _this.state.mode);
                _this.props.changeView(e, 4);
            };
            data = new FormData();
            data.append('mode', _this.state.mode);
            data.append('sequence', _this.state.sequence);
            data.append('target', _this.state.target);
            data.append('position', _this.state.position);
            data.append('replacement', _this.state.replacement);
            data.append('mutation_type', _this.state.mutation_type);
            data.append('settings', JSON.stringify(_this.formSettings));
            xhttp.send(data);
            e.preventDefault();
        };

        _this.state = {
            cursorPosition: 1,
            sequenceLength: 0,
            submitValid: false,
            loading: false,
            mode: 'PRO',
            sequence: '',
            mutation_type: '',
            target: '',
            position: 1,
            replacement: '',
            Tm_range_min: 75,
            Tm_range_max: 85,
            gc_range_min: 40,
            gc_range_max: 60,
            length_min: 25,
            length_max: 45,
            flank5_range_min: 11,
            flank5_range_max: 21,
            flank3_range_min: 11,
            flank3_range_max: 21,
            forward_overlap5: 9,
            forward_overlap3: 9,
            terminate_gc: true,
            center_mutation: true,
            primer_mode: 'complementary',
            expression_system: 'Homo sapiens'
        };
        _this.formResetDefaults = Object.assign({}, _this.state);
        return _this;
    }

    _createClass(ProteinView, [{
        key: 'render',
        value: function render() {
            var _this2 = this;

            if (this.state.loading) {
                return React.createElement(
                    'div',
                    { className: 'container mb-5 text-center' },
                    React.createElement('div', { className: 'spinner-grow mb-2', role: 'status' }),
                    React.createElement('br', null),
                    'Please wait...'
                );
            } else {
                return React.createElement(
                    'div',
                    { className: 'container mb-5' },
                    React.createElement(
                        'a',
                        {
                            className: 'btn btn-blue-grey mb-4 mr-3',
                            href: '/',
                            id: 'back',
                            onClick: function onClick(e) {
                                return _this2.props.changeView(e, 0);
                            }
                        },
                        React.createElement('i', { className: 'fas fa-arrow-left mr-3' }),
                        'main menu'
                    ),
                    React.createElement(
                        'h2',
                        { className: 'text-md-center py-2 mx-md-2 h2-responsive d-md-inline' },
                        'Protein-based Primer Design'
                    ),
                    React.createElement(
                        'form',
                        {
                            id: 'form',
                            onChange: this.formValidator,
                            onKeyUp: this.formValidator,
                            onSubmit: this.submitHandler,
                            autoComplete: 'off'
                        },
                        React.createElement(
                            'div',
                            { className: 'form-group' },
                            React.createElement(
                                'label',
                                { htmlFor: 'sequence' },
                                'Enter DNA sequence'
                            ),
                            React.createElement('i', {
                                className: 'fas fa-question-circle ml-1',
                                'data-toggle': 'tooltip',
                                title: 'Simply type in or paste your primer DNA sequence. Characters will automatically be filtered to show only A, T, C, G bases, and capitalized, as per IUPAC standards.'
                            }),
                            React.createElement('textarea', {
                                className: 'form-control md-textarea',
                                id: 'sequence',
                                name: 'sequence',
                                type: 'text',
                                value: this.state.sequence,
                                onChange: this.genericNucleobaseHandler,
                                onKeyUp: this.genericNucleobaseHandler,
                                onMouseUp: this.genericNucleobaseHandler,
                                autoFocus: true,
                                required: true,
                                style: { height: '128px' }
                            }),
                            React.createElement(
                                'div',
                                { className: 'row row-cols-1 row-cols-md-2' },
                                React.createElement(
                                    'div',
                                    { className: 'col text-left' },
                                    React.createElement(
                                        'small',
                                        { className: 'blue-grey-text' },
                                        'Cursor position: ',
                                        this.state.cursorPosition
                                    )
                                ),
                                React.createElement(
                                    'div',
                                    { className: 'col text-md-right' },
                                    React.createElement(
                                        'small',
                                        { className: 'blue-grey-text' },
                                        'Sequence length: ',
                                        this.state.sequenceLength
                                    )
                                )
                            )
                        ),
                        React.createElement(
                            'div',
                            { className: 'form-group' },
                            React.createElement(
                                'label',
                                { htmlFor: 'mutation_type' },
                                'Mutation type'
                            ),
                            React.createElement(
                                'select',
                                {
                                    className: 'browser-default custom-select',
                                    id: 'mutation_type',
                                    name: 'mutation_type',
                                    onChange: this.mutationTypeHandler,
                                    placeholder: 'Select mutation type',
                                    required: true,
                                    value: this.state.mutation_type
                                },
                                React.createElement(
                                    'option',
                                    { value: '' },
                                    'Select mutation type'
                                ),
                                React.createElement(
                                    'option',
                                    { value: 'sub' },
                                    'Substitution'
                                ),
                                React.createElement(
                                    'option',
                                    { value: 'ins' },
                                    'Insertion'
                                ),
                                React.createElement(
                                    'option',
                                    { value: 'del' },
                                    'Deletion'
                                )
                            )
                        ),
                        React.createElement(MutationCodeLayout, Object.assign({}, this.state, {
                            genericNucleobaseHandler: this.genericNucleobaseHandler,
                            genericIntHandler: this.genericIntHandler
                        })),
                        React.createElement(
                            'div',
                            { className: 'accordion my-3', id: 'accordionAdvanced' },
                            React.createElement(
                                'div',
                                { className: 'card z-depth-0 bordered' },
                                React.createElement(
                                    'div',
                                    { className: 'card-header', id: 'headingAdvanced' },
                                    React.createElement(
                                        'h5',
                                        { className: 'mb-0 text-center' },
                                        React.createElement(
                                            'button',
                                            { className: 'btn btn-link', id: 'show-advanced', type: 'button', 'data-toggle': 'collapse', 'data-target': '#advanced', 'aria-expanded': 'true', 'aria-controls': 'advanced' },
                                            'Advanced settings'
                                        )
                                    )
                                ),
                                React.createElement(
                                    'div',
                                    { id: 'advanced', className: 'collapse', 'aria-labelledby': 'headingAdvanced', 'data-parent': '#accordionAdvanced' },
                                    React.createElement(
                                        'div',
                                        { className: 'card-body' },
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3 text-center' },
                                            React.createElement('div', { className: 'col' }),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                'Min'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                'Max'
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                'Melting point (',
                                                React.createElement(
                                                    'sup',
                                                    null,
                                                    'o'
                                                ),
                                                'C)'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'Tm_range_min',
                                                        name: 'Tm_range_min',
                                                        className: 'form-control',
                                                        value: this.state.Tm_range_min,
                                                        onChange: this.genericFloatHandler
                                                    })
                                                )
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'Tm_range_max',
                                                        name: 'Tm_range_max',
                                                        className: 'form-control',
                                                        value: this.state.Tm_range_max,
                                                        onChange: this.genericFloatHandler
                                                    })
                                                )
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                'GC Content (%)'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'gc_range_min',
                                                        name: 'gc_range_min',
                                                        className: 'form-control',
                                                        value: this.state.gc_range_min,
                                                        onChange: this.genericFloatHandler
                                                    })
                                                )
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'gc_range_max',
                                                        name: 'gc_range_max',
                                                        className: 'form-control',
                                                        value: this.state.gc_range_max,
                                                        onChange: this.genericFloatHandler
                                                    })
                                                )
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                'Length (bp)'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'length_min',
                                                        name: 'length_min',
                                                        className: 'form-control',
                                                        value: this.state.length_min,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'length_max',
                                                        name: 'length_max',
                                                        className: 'form-control',
                                                        value: this.state.length_max,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                '5\' flanking region (bp)'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'flank5_range_min',
                                                        name: 'flank5_range_min',
                                                        className: 'form-control',
                                                        value: this.state.flank5_range_min,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'flank5_range_max',
                                                        name: 'flank5_range_max',
                                                        className: 'form-control',
                                                        value: this.state.flank5_range_max,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                '3\' flanking region (bp)'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'flank3_range_min',
                                                        name: 'flank3_range_min',
                                                        className: 'form-control',
                                                        value: this.state.flank3_range_min,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'flank3_range_max',
                                                        name: 'flank3_range_max',
                                                        className: 'form-control',
                                                        value: this.state.flank3_range_max,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                '5\' forward overlap (bp)'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'forward_overlap5',
                                                        name: 'forward_overlap5',
                                                        className: 'form-control',
                                                        value: this.state.forward_overlap5,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            ),
                                            React.createElement('div', { className: 'col' })
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                '3\' forward overlap (bp)'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'div',
                                                    { className: 'md-outline' },
                                                    React.createElement('input', {
                                                        type: 'number',
                                                        id: 'forward_overlap3',
                                                        name: 'forward_overlap3',
                                                        className: 'form-control',
                                                        value: this.state.forward_overlap3,
                                                        onChange: this.genericIntHandler
                                                    })
                                                )
                                            ),
                                            React.createElement('div', { className: 'col' })
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                'Expression system'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'select',
                                                    {
                                                        className: 'browser-default custom-select',
                                                        id: 'expression_system',
                                                        name: 'expression_system',
                                                        value: this.state.expression_system,
                                                        onChange: this.genericSelectHandler
                                                    },
                                                    this.props.expressionList.map(function (exp, i) {
                                                        return React.createElement(
                                                            'option',
                                                            { value: exp },
                                                            exp
                                                        );
                                                    })
                                                )
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'row row-cols-3' },
                                            React.createElement(
                                                'div',
                                                { className: 'col text-right' },
                                                'Primer type'
                                            ),
                                            React.createElement(
                                                'div',
                                                { className: 'col' },
                                                React.createElement(
                                                    'select',
                                                    {
                                                        className: 'browser-default custom-select',
                                                        id: 'primer_mode',
                                                        name: 'primer_mode',
                                                        value: this.state.primer_mode,
                                                        onChange: this.genericSelectHandler
                                                    },
                                                    React.createElement(
                                                        'option',
                                                        { value: 'complementary' },
                                                        'Complementary'
                                                    ),
                                                    React.createElement(
                                                        'option',
                                                        { value: 'overlapping' },
                                                        'Overlapping'
                                                    )
                                                )
                                            ),
                                            React.createElement('div', { className: 'col' })
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'custom-control custom-checkbox text-left mt-3' },
                                            React.createElement('input', {
                                                type: 'checkbox',
                                                className: 'custom-control-input',
                                                id: 'terminate_gc',
                                                name: 'terminate_gc',
                                                checked: this.state.terminate_gc,
                                                onChange: this.genericCheckedHandler
                                            }),
                                            React.createElement(
                                                'label',
                                                { className: 'custom-control-label', htmlFor: 'terminate_gc' },
                                                'Terminates in G/C'
                                            )
                                        ),
                                        React.createElement(
                                            'div',
                                            { className: 'custom-control custom-checkbox text-left' },
                                            React.createElement('input', {
                                                type: 'checkbox',
                                                className: 'custom-control-input',
                                                id: 'center_mutation',
                                                name: 'center_mutation',
                                                checked: this.state.center_mutation,
                                                onChange: this.genericCheckedHandler
                                            }),
                                            React.createElement(
                                                'label',
                                                { className: 'custom-control-label', htmlFor: 'center_mutation' },
                                                'Mutation at center of primer'
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        React.createElement(
                            'div',
                            { className: 'text-center mt-3' },
                            React.createElement('input', {
                                type: 'reset',
                                id: 'reset',
                                onClick: this.resetForm,
                                className: 'btn btn-warning text-dark',
                                value: 'Reset'
                            }),
                            React.createElement('input', {
                                type: 'submit',
                                className: 'btn btn-primary',
                                value: 'Submit',
                                id: 'submit',
                                disabled: !this.state.submitValid
                            })
                        )
                    )
                );
            }
        }
    }]);

    return ProteinView;
}(React.Component);