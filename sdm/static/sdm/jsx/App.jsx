const { CSSTransition } = ReactTransitionGroup;


class App extends React.Component {
    constructor(props) {
        super(props);
        this.pageId = 0;
        this.state = {
            res: null,
            pageId: 0,
            mode: null,
            transitionSpeed: 300,
            transitionName: 'fade',
        };
    }

    changeView = (e, pageId) => {
        e.preventDefault();
        this.setState({ pageId });
    }

    responseCatcher = (e, res, mode) => {
        this.setState({ res: res, mode: mode });
    }

    render() {
        if (this.state.pageId === 0) {
            return (
                <div>
                    <CSSTransition
                        key={0}
                        in={this.state.pageId === 0}
                        appear={true}
                        timeout={this.state.transitionSpeed}
                        classNames={this.state.transitionName}
                    >
                        <Menu { ...this.props } changeView={this.changeView} responseCatcher={this.responseCatcher} />
                    </CSSTransition>
                </div>
            );
        } else if (this.state.pageId === 1) {
            return (
                <div>
                    <CSSTransition
                        key={1}
                        in={this.state.pageId === 1}
                        appear={true}
                        timeout={this.state.transitionSpeed}
                        classNames={this.state.transitionName}
                    >
                        <CharacterizeView { ...this.props } changeView={this.changeView} responseCatcher={this.responseCatcher} />
                    </CSSTransition>
                </div>
            );
        } else if (this.state.pageId === 2) {
            return (
                <div>
                    <CSSTransition
                        key={2}
                        in={this.state.pageId === 2}
                        appear={true}
                        timeout={this.state.transitionSpeed}
                        classNames={this.state.transitionName}
                    >
                        <DnaView { ...this.props } changeView={this.changeView} responseCatcher={this.responseCatcher} />
                    </CSSTransition>
                </div>
            );
        } else if (this.state.pageId === 3) {
            return (
                <div>
                    <CSSTransition
                        key={3}
                        in={this.state.pageId === 3}
                        appear={true}
                        timeout={this.state.transitionSpeed}
                        classNames={this.state.transitionName}
                    >
                        <ProteinView { ...this.props } changeView={this.changeView} responseCatcher={this.responseCatcher} />
                    </CSSTransition>
                </div>
            );
        } else if (this.state.pageId === 4) {
            return (
                <div>
                    <CSSTransition
                        key={4}
                        in={this.state.pageId === 4}
                        appear={true}
                        timeout={this.state.transitionSpeed}
                        classNames={this.state.transitionName}
                    >
                        <ResultView { ...this.props } changeView={this.changeView} responseCatcher={this.responseCatcher} results={this.state.res} mode={this.state.mode} />
                    </CSSTransition>
                </div>
            );
        }
    }
}
