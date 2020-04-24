var _typeof = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function (obj) { return typeof obj; } : function (obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; };

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
            var _this2 = this;

            if (_typeof(this.props.results) === 'object') {
                if (Object.keys(this.props.results).length === 1) {
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
                } else {
                    return function () {
                        table_content = [React.createElement(
                            React.Fragment,
                            null,
                            React.createElement(
                                'button',
                                {
                                    className: 'btn btn-blue-grey mb-5',
                                    onClick: function (e) {
                                        this.props.responseCatcher(e, null);this.props.changeView(e, 0);
                                    }.bind(_this2)
                                },
                                React.createElement('i', { className: 'fas fa-arrow-left mr-2' }),
                                'main menu'
                            ),
                            React.createElement(
                                'h2',
                                { className: 'text-md-center py-2 mx-md-2 h2-responsive d-md-inline' },
                                Object.keys(_this2.props.results).length,
                                ' results'
                            )
                        )];

                        var _loop = function _loop(i, n) {
                            table_content.push(React.createElement(
                                React.Fragment,
                                null,
                                React.createElement(
                                    'div',
                                    { className: 'table-responsive text-nowrap' },
                                    React.createElement(
                                        'table',
                                        { className: 'table table-sm table-bordered table-hover' },
                                        React.createElement(
                                            'thead',
                                            null,
                                            React.createElement(
                                                'tr',
                                                null,
                                                React.createElement('th', { scope: 'col' }),
                                                React.createElement(
                                                    'th',
                                                    { scope: 'col' },
                                                    'Primer ' + (i + 1)
                                                )
                                            )
                                        ),
                                        React.createElement(
                                            'tbody',
                                            null,
                                            function () {
                                                row_content = [];
                                                for (var _key in _this2.props.results[(i + 1).toString()]) {
                                                    row_content.push(React.createElement(
                                                        React.Fragment,
                                                        null,
                                                        React.createElement(
                                                            'tr',
                                                            null,
                                                            React.createElement(
                                                                'th',
                                                                { scope: 'row' },
                                                                _key
                                                            ),
                                                            React.createElement(
                                                                'td',
                                                                null,
                                                                _this2.props.results[(i + 1).toString()][_key]
                                                            )
                                                        )
                                                    ));
                                                }
                                                return row_content;
                                            }()
                                        )
                                    )
                                )
                            ));
                        };

                        for (var i = 0, n = Object.keys(_this2.props.results).length; i < n; i++) {
                            _loop(i, n);
                        }
                        return table_content;
                    }();
                }
            } else {
                return React.createElement(
                    'div',
                    { className: 'text-center' },
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
                        'p',
                        null,
                        'An error occurred. Please try again later, or ',
                        React.createElement(
                            'a',
                            { href: 'https://github.com/kvdomingo/primerdriver/issues', target: '_blank' },
                            'report an issue'
                        ),
                        '.'
                    )
                );
            }
        }
    }]);

    return ResultView;
}(React.Component);