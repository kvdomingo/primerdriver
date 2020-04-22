class App extends React.Component {
    constructor(props) {
        super(props)
        this.state = {
            currentPage: 'main',
        }
    }

    render() {
        return (
            <div>
                <Menu {...this.props} currentPage={this.state.currentPage} />
                <CharacterizeForm {...this.props} currentPage={this.state.currentPage} />
            </div>
        );
    }
}
