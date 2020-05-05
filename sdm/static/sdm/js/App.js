var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var _ReactTransitionGroup = ReactTransitionGroup,
    CSSTransition = _ReactTransitionGroup.CSSTransition;

var App = function (_React$Component) {
    _inherits(App, _React$Component);

    function App(props) {
        _classCallCheck(this, App);

        var _this = _possibleConstructorReturn(this, (App.__proto__ || Object.getPrototypeOf(App)).call(this, props));

        _this.changeView = function (e, pageId) {
            e.preventDefault();
            _this.setState({ pageId: pageId });
        };

        _this.responseCatcher = function (e, res, mode) {
            _this.setState({ res: res, mode: mode });
        };

        _this.pageId = 0;
        _this.state = {
            res: null,
            pageId: 0,
            mode: null,
            transitionSpeed: 300,
            transitionName: 'fade'
        };
        return _this;
    }

    _createClass(App, [{
        key: 'render',
        value: function render() {
            if (this.state.pageId === 0) {
                return React.createElement(
                    'div',
                    null,
                    React.createElement(
                        CSSTransition,
                        {
                            key: 0,
                            'in': this.state.pageId === 0,
                            appear: true,
                            timeout: this.state.transitionSpeed,
                            classNames: this.state.transitionName
                        },
                        React.createElement(Menu, Object.assign({}, this.props, { changeView: this.changeView, responseCatcher: this.responseCatcher }))
                    )
                );
            } else if (this.state.pageId === 1) {
                return React.createElement(
                    'div',
                    null,
                    React.createElement(
                        CSSTransition,
                        {
                            key: 1,
                            'in': this.state.pageId === 1,
                            appear: true,
                            timeout: this.state.transitionSpeed,
                            classNames: this.state.transitionName
                        },
                        React.createElement(CharacterizeView, Object.assign({}, this.props, { changeView: this.changeView, responseCatcher: this.responseCatcher }))
                    )
                );
            } else if (this.state.pageId === 2) {
                return React.createElement(
                    'div',
                    null,
                    React.createElement(
                        CSSTransition,
                        {
                            key: 2,
                            'in': this.state.pageId === 2,
                            appear: true,
                            timeout: this.state.transitionSpeed,
                            classNames: this.state.transitionName
                        },
                        React.createElement(DnaView, Object.assign({}, this.props, { changeView: this.changeView, responseCatcher: this.responseCatcher }))
                    )
                );
            } else if (this.state.pageId === 3) {
                return React.createElement(
                    'div',
                    null,
                    React.createElement(
                        CSSTransition,
                        {
                            key: 3,
                            'in': this.state.pageId === 3,
                            appear: true,
                            timeout: this.state.transitionSpeed,
                            classNames: this.state.transitionName
                        },
                        React.createElement(ProteinView, Object.assign({}, this.props, { changeView: this.changeView, responseCatcher: this.responseCatcher }))
                    )
                );
            } else if (this.state.pageId === 4) {
                return React.createElement(
                    'div',
                    null,
                    React.createElement(
                        CSSTransition,
                        {
                            key: 4,
                            'in': this.state.pageId === 4,
                            appear: true,
                            timeout: this.state.transitionSpeed,
                            classNames: this.state.transitionName
                        },
                        React.createElement(ResultView, Object.assign({}, this.props, { changeView: this.changeView, responseCatcher: this.responseCatcher, results: this.state.res, mode: this.state.mode }))
                    )
                );
            }
        }
    }]);

    return App;
}(React.Component);