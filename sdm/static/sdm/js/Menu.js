var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var Menu = function (_React$Component) {
    _inherits(Menu, _React$Component);

    function Menu(props) {
        _classCallCheck(this, Menu);

        var _this = _possibleConstructorReturn(this, (Menu.__proto__ || Object.getPrototypeOf(Menu)).call(this, props));

        _this.pageId = 0;
        return _this;
    }

    _createClass(Menu, [{
        key: 'render',
        value: function render() {
            var _this2 = this;

            return React.createElement(
                'div',
                { className: 'row row-cols-1 row-cols-md-3' },
                this.props.stations.map(function (station, i) {
                    return React.createElement(
                        React.Fragment,
                        { key: i },
                        React.createElement(
                            'div',
                            { className: 'col' },
                            React.createElement(
                                'div',
                                { className: 'card mb-4' },
                                React.createElement(
                                    'div',
                                    { className: 'view overlay' },
                                    React.createElement('img', { className: 'card-img-top', src: station.src }),
                                    React.createElement(
                                        'a',
                                        {
                                            href: station.href,
                                            onClick: function onClick(e) {
                                                return _this2.props.changeView(e, i + 1);
                                            }
                                        },
                                        React.createElement('div', { className: 'mask rgba-black-slight' })
                                    )
                                ),
                                React.createElement(
                                    'div',
                                    { className: 'card-body text-center' },
                                    React.createElement(
                                        'a',
                                        {
                                            href: station.href,
                                            className: 'btn btn-' + station.color + ' btn-md',
                                            onClick: function onClick(e) {
                                                return _this2.props.changeView(e, i + 1);
                                            }
                                        },
                                        station.name
                                    )
                                )
                            )
                        )
                    );
                })
            );
        }
    }]);

    return Menu;
}(React.Component);