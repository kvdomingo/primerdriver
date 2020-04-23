var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var ResultView = function (_React$Component) {
    _inherits(ResultView, _React$Component);

    function ResultView(props) {
        _classCallCheck(this, ResultView);

        return _possibleConstructorReturn(this, (ResultView.__proto__ || Object.getPrototypeOf(ResultView)).call(this, props));
    }

    _createClass(ResultView, [{
        key: 'render',
        value: function render() {
            table_content = [];
            for (key in this.props.results) {
                table_content.push(React.createElement(
                    React.Fragment,
                    null,
                    React.createElement(
                        'tr',
                        null,
                        React.createElement(
                            'th',
                            { scope: 'row' },
                            key
                        ),
                        React.createElement(
                            'td',
                            null,
                            this.props.results[key]['1']
                        )
                    )
                ));
            }
            return React.createElement(
                'div',
                { className: 'table-responsive text-nowrap' },
                React.createElement(
                    'button',
                    {
                        className: 'btn btn-blue-grey mb-5',
                        onClick: function (e) {
                            this.props.responseCatcher(e, null);this.props.changeView(e, 0);
                        }.bind(this)
                    },
                    React.createElement('i', { className: 'fas fa-arrow-left mr-2' }),
                    'main menu'
                ),
                React.createElement(
                    'table',
                    { className: 'table table-sm table-bordered table-hover' },
                    React.createElement(
                        'tbody',
                        null,
                        table_content
                    )
                )
            );
        }
    }]);

    return ResultView;
}(React.Component);