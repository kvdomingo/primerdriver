class Main extends React.Component {
    constructor(props) {
        super(props)
    }

    render() {
        return (
            <div className='jumbotron blue-grey lighten-5 my-5'>
                <div className='card-deck'>
                    {this.props.stations.map((station, i) =>
                        <React.Fragment>
                            <div className='card mb-4'>
                                <div className='view overlay'>
                                    <img className='card-img-top' src={station.src} />
                                    <a href={station.href}>
                                        <div className='mask rgba-white-slight'></div>
                                    </a>
                                </div>
                                <div className='card-body text-center'>
                                    <a href={station.href} className='btn btn-md'>{station.name}</a>
                                </div>
                            </div>
                        </React.Fragment>
                    )}
                </div>
            </div>
        );
    }
}


document.addEventListener('DOMContentLoaded', () => {
    src = document.getElementsByClassName('app-data')
});
