var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var MutationCodeLayout = function (_React$Component) {
    _inherits(MutationCodeLayout, _React$Component);

    function MutationCodeLayout(props) {
        _classCallCheck(this, MutationCodeLayout);

        return _possibleConstructorReturn(this, (MutationCodeLayout.__proto__ || Object.getPrototypeOf(MutationCodeLayout)).call(this, props));
    }

    _createClass(MutationCodeLayout, [{
        key: 'render',
        value: function render() {
            if (this.props.mutation_type === 'sub') {
                return React.createElement(
                    'div',
                    { className: 'row row-cols-1 row-cols-md-3' },
                    React.createElement(
                        'div',
                        { className: 'col form-group' },
                        React.createElement(
                            'label',
                            { htmlFor: 'target' },
                            'Target'
                        ),
                        React.createElement('input', {
                            type: 'text',
                            className: 'form-control',
                            name: 'target',
                            id: 'target',
                            value: this.props.target,
                            onChange: this.props.genericNucleobaseHandler,
                            required: true
                        })
                    ),
                    React.createElement(
                        'div',
                        { className: 'col form-group' },
                        React.createElement(
                            'label',
                            { htmlFor: 'target' },
                            'Position'
                        ),
                        React.createElement('input', {
                            type: 'number',
                            min: '0',
                            className: 'form-control',
                            name: 'position',
                            id: 'position',
                            value: this.props.position,
                            onChange: this.props.genericIntHandler,
                            required: true
                        })
                    ),
                    React.createElement(
                        'div',
                        { className: 'col form-group' },
                        React.createElement(
                            'label',
                            { htmlFor: 'target' },
                            'Replacement'
                        ),
                        React.createElement('input', {
                            type: 'text',
                            className: 'form-control',
                            name: 'replacement',
                            id: 'replacement',
                            min: '0',
                            value: this.props.replacement,
                            onChange: this.props.genericNucleobaseHandler,
                            required: true
                        })
                    )
                );
            } else if (this.props.mutation_type === 'ins') {
                return React.createElement(
                    'div',
                    { className: 'row row-cols-1 row-cols-md-2' },
                    React.createElement(
                        'div',
                        { className: 'col form-group' },
                        React.createElement(
                            'label',
                            { htmlFor: 'target' },
                            'Position'
                        ),
                        React.createElement('input', {
                            type: 'number',
                            min: '0',
                            className: 'form-control',
                            name: 'position',
                            id: 'position',
                            value: this.props.position,
                            onChange: this.props.genericIntHandler,
                            required: true
                        })
                    ),
                    React.createElement(
                        'div',
                        { className: 'col form-group' },
                        React.createElement(
                            'label',
                            { htmlFor: 'target' },
                            'Replacement'
                        ),
                        React.createElement('input', {
                            type: 'text',
                            className: 'form-control',
                            name: 'replacement',
                            id: 'replacement',
                            min: '0',
                            value: this.props.replacement,
                            onChange: this.props.genericNucleobaseHandler,
                            required: true
                        })
                    )
                );
            } else if (this.props.mutation_type === 'del') {
                return React.createElement(
                    'div',
                    { className: 'row row-cols-1 row-cols-md-2' },
                    React.createElement(
                        'div',
                        { className: 'col form-group' },
                        React.createElement(
                            'label',
                            { htmlFor: 'target' },
                            'Target'
                        ),
                        React.createElement('input', {
                            type: 'text',
                            className: 'form-control',
                            name: 'target',
                            id: 'target',
                            value: this.props.target,
                            onChange: this.props.genericNucleobaseHandler,
                            required: true
                        })
                    ),
                    React.createElement(
                        'div',
                        { className: 'col form-group' },
                        React.createElement(
                            'label',
                            { htmlFor: 'target' },
                            'Position'
                        ),
                        React.createElement('input', {
                            type: 'number',
                            min: '0',
                            className: 'form-control',
                            name: 'position',
                            id: 'position',
                            value: this.props.position,
                            onChange: this.props.genericIntHandler,
                            required: true
                        })
                    )
                );
            } else {
                return null;
            }
        }
    }]);

    return MutationCodeLayout;
}(React.Component);