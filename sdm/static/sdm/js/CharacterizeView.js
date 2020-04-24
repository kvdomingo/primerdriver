var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var CharacterizeView = function (_React$Component) {
    _inherits(CharacterizeView, _React$Component);

    function CharacterizeView(props) {
        _classCallCheck(this, CharacterizeView);

        var _this = _possibleConstructorReturn(this, (CharacterizeView.__proto__ || Object.getPrototypeOf(CharacterizeView)).call(this, props));

        _this.textAreaHandler = function (e) {
            sequence = e.target.value.toUpperCase().split("");
            filteredSequence = [];
            sequence.forEach(function (char, i) {
                if (['A', 'T', 'C', 'G'].includes(char)) filteredSequence.push(char);
            });
            _this.setState({
                sequence: filteredSequence.join(""),
                cursorPosition: e.target.selectionStart + 1,
                sequenceLength: e.target.value.length
            });
        };

        _this.mismatchedBasesHandler = function (e) {
            _this.setState({ mismatched_bases: parseInt(e.target.value) });
        };

        _this.mutationTypeHandler = function (e) {
            _this.setState({ mutation_type: e.target.value });
        };

        _this.resetForm = function () {
            _this.setState(Object.assign({}, _this.formResetDefaults));
        };

        _this.formValidator = function () {
            validSequence = _this.state.sequenceLength > 0;
            validMismatch = _this.state.mismatched_bases > 0;
            validSequence && validMismatch ? _this.setState({ submitValid: true }) : _this.setState({ submitValid: false });
            _this.formData = {
                mode: _this.state.mode,
                sequence: _this.state.sequence,
                mismatched_bases: _this.state.mismatched_bases,
                mutation_type: _this.state.mutation_type
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
                }
                _this.props.responseCatcher(e, res, _this.state.mode);
                _this.props.changeView(e, 4);
            };
            data = new FormData();
            data.append('mode', _this.state.mode);
            data.append('sequence', _this.state.sequence);
            data.append('mismatched_bases', _this.state.mismatched_bases);
            data.append('mutation_type', _this.state.mutation_type);
            xhttp.send(data);
            e.preventDefault();
        };

        _this.state = {
            cursorPosition: 1,
            sequenceLength: 0,
            submitValid: false,
            loading: false,
            mode: 'CHAR',
            sequence: '',
            mismatched_bases: 0,
            mutation_type: ''
        };
        _this.formResetDefaults = Object.assign({}, _this.state);
        return _this;
    }

    _createClass(CharacterizeView, [{
        key: 'render',
        value: function render() {
            var _this2 = this;

            if (this.state.loading) {
                return React.createElement(
                    'div',
                    { className: 'container mb-5' },
                    React.createElement(
                        'div',
                        { className: 'spinner-grow', role: 'status' },
                        React.createElement(
                            'span',
                            { className: 'sr-only' },
                            'Loading...'
                        )
                    )
                );
            } else {
                return React.createElement(
                    'div',
                    { className: 'container mb-5' },
                    React.createElement(
                        'a',
                        {
                            className: 'btn btn-blue-grey btn-rounded mb-4',
                            href: '/',
                            id: 'back',
                            onClick: function onClick(e) {
                                return _this2.props.changeView(e, 0);
                            }
                        },
                        React.createElement('i', { className: 'fas fa-arrow-left mr-2' }),
                        'main menu'
                    ),
                    React.createElement(
                        'h2',
                        { className: 'text-md-center py-2 mx-md-2 h2-responsive d-md-inline' },
                        'Characterization'
                    ),
                    React.createElement(
                        'form',
                        {
                            id: 'form',
                            onChange: this.formValidator,
                            onSubmit: this.submitHandler
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
                                onChange: this.textAreaHandler,
                                onKeyUp: this.textAreaHandler,
                                onMouseUp: this.textAreaHandler,
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
                            { className: 'row row-cols-1 row-cols-md-2' },
                            React.createElement(
                                'div',
                                { className: 'col form-group' },
                                React.createElement(
                                    'label',
                                    { htmlFor: 'mismatched_bases' },
                                    'Number of mismatched bases'
                                ),
                                React.createElement('input', {
                                    className: 'form-control',
                                    min: '0',
                                    type: 'number',
                                    id: 'mismatched_bases',
                                    name: 'mismatched_bases',
                                    value: this.state.mismatched_bases,
                                    onChange: this.mismatchedBasesHandler,
                                    required: true
                                })
                            ),
                            React.createElement(
                                'div',
                                { className: 'col form-group' },
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
                                        { value: '', disabled: true },
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

    return CharacterizeView;
}(React.Component);