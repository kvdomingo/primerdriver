var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var Main = function (_React$Component) {
    _inherits(Main, _React$Component);

    function Main(props) {
        _classCallCheck(this, Main);

        return _possibleConstructorReturn(this, (Main.__proto__ || Object.getPrototypeOf(Main)).call(this, props));
    }

    _createClass(Main, [{
        key: 'render',
        value: function render() {
            return React.createElement(
                'div',
                { className: 'jumbotron blue-grey lighten-5 my-5' },
                React.createElement(
                    'div',
                    { className: 'card-deck' },
                    this.props.stations.map(function (station, i) {
                        return React.createElement(
                            React.Fragment,
                            null,
                            React.createElement(
                                'div',
                                { className: 'card mb-4' },
                                React.createElement(
                                    'div',
                                    { className: 'view overlay' },
                                    React.createElement('img', { className: 'card-img-top', src: station.src }),
                                    React.createElement(
                                        'a',
                                        { href: station.href },
                                        React.createElement('div', { className: 'mask rgba-white-slight' })
                                    )
                                ),
                                React.createElement(
                                    'div',
                                    { className: 'card-body text-center' },
                                    React.createElement(
                                        'a',
                                        { href: station.href, className: 'btn btn-md' },
                                        station.name
                                    )
                                )
                            )
                        );
                    })
                )
            );
        }
    }]);

    return Main;
}(React.Component);

document.addEventListener('DOMContentLoaded', function () {
    src = document.getElementsByClassName('app-data');
});