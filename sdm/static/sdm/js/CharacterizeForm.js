var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var CharacterizeForm = function (_React$Component) {
    _inherits(CharacterizeForm, _React$Component);

    function CharacterizeForm(props) {
        _classCallCheck(this, CharacterizeForm);

        var _this = _possibleConstructorReturn(this, (CharacterizeForm.__proto__ || Object.getPrototypeOf(CharacterizeForm)).call(this, props));

        if (_this.props.currentPage !== 'characterize') {
            var _ret;

            return _ret = null, _possibleConstructorReturn(_this, _ret);
        }
        return _this;
    }

    _createClass(CharacterizeForm, [{
        key: 'render',
        value: function render() {
            return React.createElement(
                'div',
                { className: 'form-group' },
                React.createElement(
                    'label',
                    { htmlFor: 'sequence' },
                    'Enter DNA sequence'
                ),
                React.createElement('textarea', {
                    className: 'form-control md-textarea',
                    id: 'sequence',
                    name: 'sequence',
                    type: 'text',
                    value: this.props.sequence,
                    onChange: this.props.handleChange,
                    autoFocus: true
                })
            );
        }
    }]);

    return CharacterizeForm;
}(React.Component);